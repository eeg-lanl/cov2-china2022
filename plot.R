require(tidyverse)
require(lubridate)
require(coda)
require(patchwork)
require(RColorBrewer)

get_hpd <- function(x)
{
    x <- na.omit(x)
    ans0 <- median(x)
    if (is.na(ans0))
    {
        ans1 <- rep(NA, 2)
        ans2 <- rep(NA, 2)
    } else {
        ans1 <- HPDinterval(mcmc(x), prob=0.95)
        ans2 <- HPDinterval(mcmc(x), prob=0.50)
    }
    ans <- data.frame(mid=ans0,
                      low1=ans1[1], high1=ans1[2],
                      low2=ans2[1], high2=ans2[2])
    return(ans)
}

calculate_r_R <- function(beta)
{
    cal_R0 <- function(rr)
    {
        func_R0 <- function(tau)
            dgamma(x=tau, shape=3, scale=3.3/3, log=FALSE)*exp(-rr*tau)
        res <- integrate(func_R0, lower=0, upper=Inf)
        return(1 / res$value)
    }

    k_E <- 2 / 1.9
    k <- 1 / 1.5
    f <- 0.913765 # median est
    JacM <- matrix(c(-k_E, 0, beta, beta, beta, beta,
                     k_E, -k_E, 0, 0, 0, 0,
                     0, f*k_E, -k, 0, 0, 0,
                     0, (1-f)*k_E, 0, -k, 0, 0,
                     0, 0, k, 0, -k, 0,
                     0, 0, 0, k, 0, -k),
                     6, 6, byrow=TRUE)
    res <- eigen(JacM)
    rr <- max(Re(res$values))
    r0 <- cal_R0(rr)

    return(data.frame("r" = rr, "R" = r0))
}

#--------------------------------------------------
# Setup
#--------------------------------------------------

pc <- c("TRUE" = "black", "FALSE" = "#777777",
        "survS" = "#5c3566", "survI" = "#cc0000", "survR" = "#4e9a06",
        "S" = "#3465a4", "E" = "#c17d11", "I" = "#f57900", "R" = "#73d216")

model <- "seir"
label <- "base"
load(paste0(model, "-", label, ".rda"))

dec26 <- tibble(Name = c("survS", "survI", "survR"), Value = data_list$dec26)

#--------------------------------------------------
# Cases with model fit
#--------------------------------------------------

cases <- read_csv("cases.csv") %>%
         rename(Day = ymd, Province = province) %>%
         filter(!(Province %in% c("aomen", "taiwan", "xianggang"))) %>%
         filter(Day >= ymd("2022-10-28"), Day <= ymd("2022-12-26")) %>%
         select(c(Day, Province, locAsymNum, locIncrNum)) %>%
         mutate(Cases = locAsymNum + locIncrNum) %>%
         group_by(Day) %>%
         summarize(Cases = sum(Cases)) %>%
         ungroup() %>%
         arrange(Day) %>%
         mutate(Used = Day <= ymd("2022-11-11"))

y14 <- cbind(chain$y1[,,14], chain$y2[,,14], chain$y3[,,14])
c14 <- y14[, 2:ncol(y14)] - y14[, 1:(ncol(y14)-1)]
m_cases <- bind_rows(apply(c14, 2, get_hpd)) %>%
           as_tibble() %>%
           mutate(Day=seq(ymd("2022-10-24"), ymd("2022-12-26"), by=1))

dat_cases <- full_join(cases, m_cases) %>%
             filter(Day >= ymd("2022-10-28"))

p1 <- ggplot(data=dat_cases) +
      geom_ribbon(aes(x=Day, ymin=low1, ymax=high1), color=NA, fill="#fcaf3e", alpha=0.5) +
      geom_ribbon(aes(x=Day, ymin=low2, ymax=high2), color=NA, fill="#fcaf3e", alpha=0.75) +
      geom_line(aes(x=Day, y=mid), color="orange2", size=1.5) +
      geom_vline(data=data.frame(x=c(ymd("2022-11-11", "2022-12-07"))), aes(xintercept=x)) +
      geom_point(aes(x=Day, y=Cases, color=Used)) +
      scale_y_continuous(trans="log10") +
      scale_color_manual(values=pc) +
      guides(color="none") +
      theme_bw() +
      theme(axis.title.x=element_blank()) +
      theme(text=element_text(size=16)) +
      labs(y="daily cases")

#--------------------------------------------------
# SEIR
#--------------------------------------------------

get_seir <- function(y, day1, dayN)
{
    add_cols <- function(x0, lab)
    {
        bind_rows(apply(x0, 2, function(x) get_hpd(x))) %>%
            add_column(day=seq(ymd(day1), ymd(dayN), by=1)) %>%
            add_column(epi_state=lab)
    }
    model_S <- add_cols(y[,,1], "S")
    model_E <- add_cols(rowSums(y[,,2:3], dim=2), "E")
    model_I <- add_cols(rowSums(y[,,c(4:5, 10:11)], dim=2), "I")
    model_R <- add_cols(rowSums(y[,,c(6:9, 12)], dim=2), "R")
    model_seir <- bind_rows(model_S, model_E, model_I, model_R) %>%
                  as_tibble() %>%
                  mutate(epi_state=factor(epi_state, levels=c("S", "E", "I", "R")))
    model_seir
}

model_seir1 <- get_seir(chain$y1, "2022-10-23", "2022-11-11")
model_seir2 <- get_seir(chain$y2, "2022-11-12", "2022-12-07")
model_seir3 <- get_seir(chain$y3, "2022-12-08", "2022-12-26")
model_seir <- bind_rows(model_seir1, model_seir2, model_seir3)

p2 <- ggplot(data=model_seir, aes(x=day)) +
      geom_ribbon(aes(ymin=low1, ymax=high1, fill=epi_state), alpha=0.5) +
      geom_ribbon(aes(ymin=low2, ymax=high2, fill=epi_state), alpha=0.5) +
      geom_line(aes(y=mid, color=epi_state), size=1.2) +
      scale_color_manual(values=pc, limits=force) +
      scale_fill_manual(values=pc, limits=force) +
      labs(y="prevalence", color="", fill="") +
      theme_bw() +
      theme(text=element_text(size=16)) +
      theme(axis.title.x=element_blank())

#--------------------------------------------------
# Dec 26
#--------------------------------------------------

dat_m <- as_tibble(chain[["m_dec26"]], .name_repair=~c("survS", "survI", "survR")) %>%
         pivot_longer(everything()) %>%
         mutate(name=factor(name, levels=c("survS", "survI", "survR")))
dat_p <- as_tibble(chain[["ppd_dec26"]], .name_repair=~c("survS", "survI", "survR")) %>%
         pivot_longer(everything()) %>%
         mutate(name=factor(name, levels=c("survS", "survI", "survR")))

p3 <- ggplot() +
      geom_density(data=dat_p, aes(x=value, fill=name), color=NA, alpha=0.4, show.legend=F) +
      geom_density(data=dat_m, aes(x=value, color=name), size=1.1) +
      geom_vline(data=dec26, aes(xintercept=Value, color=Name), size=1.1, linetype="dashed", show.legend=F) +
      labs(x="prevalence", y="prob density") +
      scale_color_manual(values=pc, limits=force) +
      scale_fill_manual(values=pc, limits=force) +
      theme_bw() +
      theme(text=element_text(size=16)) +
      theme(legend.title=element_blank()) +
      theme(axis.text.y=element_blank())

#--------------------------------------------------
# Parameters
#--------------------------------------------------

model_beta <- data.frame(beta1=chain[["beta1"]], beta2=chain[["beta2"]], beta3=chain[["beta3"]]) %>%
              pivot_longer(everything())

trans_hpd <- model_beta %>% group_by(name) %>%
             summarize(hpd = get_hpd(value)) %>% unpack(hpd) %>%
             ungroup() %>% relocate(mid, .after=4) %>%
             rename(period = name) %>%
             pivot_longer(-1) %>%
             rename(beta = value, limit = name) %>%
             rowwise() %>%
             mutate(rR = calculate_r_R(beta)) %>% unpack(rR) %>%
             rename("r (exponential)" = "r", "R (exponential)" = "R") %>%
             pivot_longer(-c(1,2)) %>%
             pivot_wider(names_from = limit) %>%
             mutate(daterange = case_when(period=="beta1" ~ "Zero COVID",
                                          period=="beta2" ~ "20 Measures",
                                          period=="beta3" ~ "10 Measures"))

p4 <- ggplot(data=filter(trans_hpd, name!="beta"), aes(y=daterange)) +
      geom_point(aes(x=mid), size=3) +
      geom_linerange(aes(xmin=low1, xmax=high1), size=0.5) +
      geom_linerange(aes(xmin=low2, xmax=high2), size=1.5) +
      facet_wrap(~name, ncol=1, scales="free_x") +
      theme_bw() +
      theme(text=element_text(size=16)) +
      theme(axis.title=element_blank())

#--------------------------------------------------
# Composite
#--------------------------------------------------

p2a <- p2 + coord_cartesian(x=c(ymd("2022-12-07"), NA)) + theme(legend.position=c(0.1, 0.6))
p3a <- p3 + theme(legend.position=c(0.83, 0.6)) + coord_cartesian(x=c(2.6e8, 7.7e8))

top <- p1 + p3a + plot_layout(widths=c(2,1))
bot <- p4 + p2a + plot_layout(widths=c(1,1.8))

pdf(file="data_seir.pdf", width=9.5, height=7.5)
wrap_elements(full=top) + wrap_elements(full=bot) + plot_layout(ncol=1, heights=c(0.8, 1))
dev.off()

#--------------------------------------------------
# Compare with National survey
#--------------------------------------------------

nat1 <- read_csv("national_pcr.csv")
nat2 <- read_csv("national_ant.csv")
nat <- bind_rows(nat1, nat2) %>%
       filter(Day <= ymd("2023-01-22"))

y14 <- cbind(chain$y1[,,14], chain$y2[,,14], chain$y3[,,14], chain$y4[,,14])
c14 <- y14[, 2:ncol(y14)] - y14[, 1:(ncol(y14)-1)]

yR <- cbind(chain$y1[,,9], chain$y2[,,9], chain$y3[,,9], chain$y4[,,9]) +
        cbind(chain$y1[,,12], chain$y2[,,12], chain$y3[,,12], chain$y4[,,12])
R <- data.frame(Day=seq(ymd("2022-10-23"), ymd("2023-01-22"), by=1), R=apply(yR, 2, median))

m_cases <- bind_rows(apply(c14, 2, get_hpd)) %>%
           as_tibble() %>%
           mutate(Day=seq(ymd("2022-10-24"), ymd("2023-01-22"), by=1)) %>%
           left_join(R) %>%
           filter(Day >= min(nat$Day), Day <= max(nat$Day))

dat_nat <- nat %>%
           left_join(m_cases) %>%
           mutate(w = case_when(Type == "PCR" ~ Tests / ((data_list$N-R) * 0.444),
                                Type == "ANT" ~ Tests / ((data_list$N-R) * 0.522))) %>%
           mutate(across(mid:high2, ~ .x * w))

p <- ggplot(data=m_cases, aes(x=Day)) +
     geom_vline(aes(xintercept=ymd("2022-12-26")), linetype="dashed") +
     # model prediction of full cases:
     geom_ribbon(aes(ymin=low1, ymax=high1), alpha=0.5) +
     geom_ribbon(aes(ymin=low2, ymax=high2), alpha=0.75) +
     geom_line(aes(y=mid), size=1.3) +
     # model prediction with testing correction:
     geom_ribbon(data=dat_nat, aes(ymin=low1, ymax=high1, fill=Type), alpha=0.5) +
     geom_ribbon(data=dat_nat, aes(ymin=low2, ymax=high2, fill=Type), alpha=0.75) +
     geom_line(data=dat_nat, aes(y=mid, color=Type), size=1.3) +
     # data from PCR and ANT tests:
     geom_point(data=dat_nat, aes(y=Cases, fill=Type), pch=21, size=2, stroke=1.5) +
     scale_x_date(limits=c(min(dat_nat$Day), ymd("2022-12-31"))) +
     scale_y_continuous(trans="log10", limits=c(1e3, NA)) +
     labs(y = "cases") +
     theme_bw() +
     theme(axis.title.x=element_blank()) +
     theme(text=element_text(size=14)) +
     theme(legend.position=c(0, 1), legend.justification=c(0, 1))

pdf(file="nat_test.pdf", width=5, height=4)
p
dev.off()

### least-squared fit to get the "voluntary testing fraction" used above ###

#dat <- nat %>% filter(Type == "PCR") %>%
dat <- nat %>% filter(Type == "ANT") %>%
       left_join(m_cases) %>%
       filter(Day <= ymd("2022-12-31"))
myfunc <- function(v)
    sum((dat$mid * dat$Tests / ((data_list$N-dat$R) * v) - dat$Cases)^2)
optimize(myfunc, interval=c(0.01, 1))
