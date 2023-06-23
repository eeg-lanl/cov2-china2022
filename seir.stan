functions {
    vector seir(real t, vector y, real beta, real f, int N) {

        real S = y[1];
        real I = y[4] + y[5] + y[10] + y[11];

        vector[14] dy_dt;

        real kE = 2 / 1.9;
        real kA = 1 / 1.5;
        real kIP = 1 / 1.5;
        real kIS = 1 / 1.5;
        real kS = 3 / 4.2;
        real w = 1;
        real kT = 0.5;

        dy_dt[1] = -beta * I * S / N;
        dy_dt[2] =  beta * I * S / N - kE * y[2];
        dy_dt[3] = kE * y[2] - kE * y[3];
        dy_dt[4] = kE * f * y[3] - kIP * y[4];
        dy_dt[5] = kIP * y[4] - kIS * y[5];
        dy_dt[6] = kIS * y[5] - kS * y[6];
        dy_dt[7] = kS * y[6] - kS * y[7];
        dy_dt[8] = kS * y[7] - kS * y[8];
        dy_dt[9] = kS * y[8];
        dy_dt[10] = kE * (1-f) * y[3] - kA * y[10];
        dy_dt[11] = kA * y[10] - kA * y[11];
        dy_dt[12] = kA * y[11];
        dy_dt[13] = kE * w * y[3] - kT * y[13];
        dy_dt[14] = kT * y[13];

        return dy_dt;
    }
}

data {
    int<lower=0> N;                  // total population size
    array[5] int n_days;             // length of each ODE period, length of case data period
    array[n_days[1]] real days1;     // beta1, Oct 22 - Nov 11
    array[n_days[2]] real days2;     // beta2, Nov 11 - Dec 07
    array[n_days[3]] real days3;     // beta3, Dec 07 - Dec 26
    array[n_days[4]] int daysC;      // days with case data
    array[n_days[4]] int cases;      // daily cases, Oct 28 - Nov 11
    array[n_days[5]] real days4;     // beta3, Dec 26 - Jan 22
    array[3] int dec26;              // survey data
}

parameters {
    real<lower=1> E0; // for the initial ODE states
    real<lower=0, upper=1> f;
    real<lower=0> phi_inv;
    real<lower=0> beta1;
    real<lower=0> beta3;
    real<lower=0> beta_diff21; // to enforce beta2 > beta1
}

transformed parameters {
    real phi = 1. / phi_inv;
    real beta2 = beta1 + beta_diff21;
    vector[14] init;

    /* integrate for each time period */

    init = [N-(E0*2), E0, E0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    array[n_days[1]-1] vector<lower=0>[14] y1 = ode_rk45(seir, init, days1[1], days1[2:], beta1, f, N);
    array[n_days[2]-1] vector<lower=0>[14] y2 = ode_rk45(seir, y1[size(y1)], days2[1], days2[2:], beta2, f, N);
    array[n_days[3]-1] vector<lower=0>[14] y3 = ode_rk45(seir, y2[size(y2)], days3[1], days3[2:], beta3, f, N);
    array[n_days[5]-1] vector<lower=0>[14] y4 = ode_rk45(seir, y3[size(y3)], days4[1], days4[2:], beta3, f, N);

    /* model-expected value for each day of data */

    vector[n_days[4]] m_cases;
    for (i in 1:n_days[4])
    {
        m_cases[i] = y1[daysC[i]-1, 14] - y1[daysC[i]-2, 14];
    }

    vector[3] m_dec26;
    m_dec26[1] = sum(y3[size(y3)][1:4]) + sum(y3[size(y3)][10:12]); // report never been infected
    m_dec26[2] = sum(y3[size(y3)][5:8]);                            // report infected, because they have symptoms
    m_dec26[3] = y3[size(y3)][9];                                   // report recovered, because they used to have symptoms
}

model {
    /* priors */
    E0 ~ uniform(1, cases[1]*10);
    f ~ uniform(0, 1);
    beta1 ~ normal(3, 2);
    beta3 ~ normal(3, 2);
    beta_diff21 ~ exponential(1);
    phi_inv ~ exponential(5);

    /* likelihood, sampling distribution */
    cases ~ neg_binomial_2(m_cases, phi);
    for (i in 1:3) {
        dec26[i] ~ neg_binomial_2(m_dec26[i], phi);
    }
}

generated quantities {
    array[n_days[4]] int ppd_cases = neg_binomial_2_rng(m_cases, phi);
    vector[3] ppd_dec26;
    for (i in 1:3) {
        ppd_dec26[i] = neg_binomial_2_rng(m_dec26[i], phi);
    }
}
