####load library####
library(R2jags)
library(tidyverse)
library(dplyr)

####load Simulations_method with_GS_interaction for all functions ####

#### function for 2 subset for HIgh data with Millen method######
####functions#####
stopping_rules_first_analyse_variation_2_subset_millen <-
  function(P1_0, P1_1, P5_0, P5_1, gamma_1, epsilon_1) {
    nDecision <- "continue with 2 subsets"
    efficacy_1 <- ifelse(P1_1 > gamma_1, 1, 0)
    efficacy_0 <- ifelse(P1_0 > gamma_1, 1, 0)
    interaction_1 <- ifelse(P5_1 > epsilon_1, 1, 0)
    interaction_0 <- ifelse(P5_0 > epsilon_1, 1, 0)
    enrichment_1 <-
      ifelse(efficacy_1 == 1 & interaction_0 == 1, 1, 0)
    enrichment_0 <-
      ifelse(efficacy_0 == 1 & interaction_1 == 1, 1, 0)
    continue <- 0
    if (efficacy_1 == 1 & interaction_0 == 1) {
      nDecision <- "efficacy in the subset 1"
    } else if (efficacy_0 == 1 & interaction_1 == 1) {
      nDecision <- "efficacy in the subset 0"
    } else{
      nDecision <- "continue with 2 subsets"
      continue <- 1
    }
    return(
      c(
        nDecision,
        continue,
        efficacy_0,
        efficacy_1,
        interaction_0,
        interaction_1,
        enrichment_0,
        enrichment_1
      )
    )
  }

stopping_rules_other_analyse_variation_2_subset_millen <-
  function(nDecision,
           P1_0,
           P1_1,
           P5_0,
           P5_1,
           gamma_1,
           epsilon_1) {
    continue <- 0
    efficacy_1 <- 0
    efficacy_0 <- 0
    interaction_1 <- 0
    interaction_0 <- 0
    enrichment_1 <- 0
    enrichment_0 <- 0
    if (nDecision == "efficacy in the subset 1") {
      nDecision <- "efficacy in the subset 1"
      enrichment_1 <- 1
      efficacy_1 <- 1
      interaction_0 <- 1
    }
    else if (nDecision == "efficacy in the subset 0") {
      nDecision <- "efficacy in the subset 0"
      enrichment_0 <- 1
      efficacy_0 <- 1
      interaction_1 <- 1
    }
    else {
      efficacy_1 <- ifelse(P1_1 > gamma_1, 1, 0)
      efficacy_0 <- ifelse(P1_0 > gamma_1, 1, 0)
      interaction_1 <- ifelse(P5_1 > epsilon_1, 1, 0)
      interaction_0 <- ifelse(P5_0 > epsilon_1, 1, 0)
      enrichment_1 <-
        ifelse(efficacy_1 == 1 & interaction_0 == 1, 1, 0)
      enrichment_0 <-
        ifelse(efficacy_0 == 1 & interaction_1 == 1, 1, 0)
      if (efficacy_1 == 1 & interaction_0 == 1) {
        nDecision <- "efficacy in the subset 1"
      }
      else if (efficacy_0 == 1 & interaction_1 == 1) {
        nDecision <- "efficacy in the subset 0"
      } else{
        nDecision <- "continue with 2 subsets"
        continue <- 1
      }
    }
    return(
      c(
        nDecision,
        continue,
        efficacy_0,
        efficacy_1,
        interaction_0,
        interaction_1,
        enrichment_0,
        enrichment_1
      )
    )
  }


create_matrix_HIGH <- function() {
  mat.sub1 <-
    matrix(nrow = 1,
           ncol = 45,
           dimnames = list(
             NULL,
             c(
               "RR.global",
               "2.5%",
               "97.5%",
               "RR.subset1",
               "2.5%",
               "97.5%",
               "RR.subset0",
               "2.5%",
               "97.5%",
               "continue_1 with subset ",
               "efficacy0_1",
               "efficacy1_1",
               "interaction0_1",
               "interaction1_1",
               "enrichment0_1",
               "enrichment1_1",
               "continue_2 with subset ",
               "efficacy0_2",
               "efficacy1_2",
               "interaction0_2",
               "interaction1_2",
               "enrichment0_2",
               "enrichment1_2",
               "continue_3 with subset ",
               "efficacy0_3",
               "efficacy1_3",
               "interaction0_3",
               "interaction1_3",
               "enrichment0_3",
               "enrichment1_3",
               "continue_4 with subset ",
               "efficacy0_4",
               "efficacy1_4",
               "interaction0_4",
               "interaction1_4",
               "enrichment0_4",
               "enrichment1_4",
               "N_subset1_stage1",
               "N_subset0_stage1",
               "N_subset1_stage2",
               "N_subset0_stage2",
               "N_subset1_stage3",
               "N_subset0_stage3",
               "N_subset1_stage4",
               "N_subset0_stage4"
             )
           ))
}# create matrix for simulation


binary_bayesian_model_2_subset_Millen <- function() {
  nb11 ~ dbin(p11, n11)
  nb10 ~ dbin(p10, n10)
  nb01 ~ dbin(p01, n01)
  nb00 ~ dbin(p00, n00)
  
  p10 ~ dbeta(1, 1)
  p00 ~ dbeta(1, 1)
  p01 ~ dbeta(1, 1)
  p11 ~ dbeta(1, 1)
  
  
  p1.estim <- (p11 * pa) + (p10 * (1 - pa)) ### probability in treatment group
  p0.estim <- (p01 * pa) + (p00 * (1 - pa))  ### probability in control group
  RR <- p1.estim / p0.estim
  RR.subset1 <- p11 / p01 # effect of treatment in subset 1
  RR.subset0 <- p10 / p00 # effect of treatment in subset 0
  interaction_0 <- RR.subset0 / RR.subset1
  interaction_1 <- RR.subset1 / RR.subset0
  P1_1 <- (RR.subset1 < lambda_1)
  P1_0 <- (RR.subset0 < lambda_1)
  P5_0 <- (interaction_0 > eta_1)
  P5_1 <- (interaction_1 > eta_1)
  biais_0 <- RR.subset0 - (pi10 / pi00)
  biais_1 <- RR.subset1 - (pi11 / pi01)
  biais_int_0 <-
    interaction_0 - ((pi10 / pi00) / (pi11 / pi01))
  biais_int_1 <-
    interaction_1 - ((pi11 / pi01) / (pi10 / pi00))
  
} # bayesian model set 1 for jags function


create_list_for_inits_values <-
  function(Data,
           pi11,
           pi10,
           pi01,
           pi00,
           pa,
           lambda_1,
           eta_1) {
    list(
      n11 = sum(Data$group == 1 & Data$subset == 1),
      n10 = sum(Data$group == 1 &
                  Data$subset == 0),
      n01 = sum(Data$group == 0 &
                  Data$subset == 1),
      n00 = sum(Data$group == 0 &
                  Data$subset == 0),
      nb11 = sum(Data$group == 1 & Data$subset == 1 & Data$y == 1,
                 na.rm = T),
      nb10 = sum(Data$group == 1 & Data$subset == 0 & Data$y == 1,
                 na.rm = T),
      nb01 = sum(Data$group == 0 & Data$subset == 1 & Data$y == 1,
                 na.rm = T),
      nb00 = sum(Data$group == 0 & Data$subset == 0 & Data$y == 1,
                 na.rm = T),
      pi11 = pi11,
      pi10 = pi10,
      pi01 = pi01,
      pi00 = pi00,
      pa = pa,
      lambda_1 = lambda_1,
      eta_1 = eta_1
    )
  }# list for iniatialization values for jags

results_bayesian_model_HIGH_DATA <-
  function(bayesian_jags_fonction_IA1, Data) {
    p11 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p11
    p10 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p10
    p01 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p01
    p00 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p00
    RR <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR
    RR_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[5, 3]
    RR_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[5, 7]
    RR.subset1 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR.subset1
    RR.subset1_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[7, 3]
    RR.subset1_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[7, 7]
    RR.subset0 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR.subset0
    RR.subset0_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[6, 3]
    RR.subset0_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[6, 7]
    P1_0 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P1_0
    P1_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P1_1
    P5_0 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P5_0
    P5_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P5_1
    interaction_0 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$interaction_0
    interaction_1 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$interaction_1
    N_subset1_stage1 <- sum(Data$subset == 1)
    N_subset0_stage1 <- sum(Data$subset == 0)
    
    results <-
      cbind(
        RR,
        RR_2.5,
        RR_97.5,
        RR.subset1,
        RR.subset1_2.5,
        RR.subset1_97.5,
        RR.subset0,
        RR.subset0_2.5,
        RR.subset0_97.5,
        P1_0,
        P1_1,
        P5_0,
        P5_1,
        N_subset1_stage1,
        N_subset0_stage1
      )
    results
  }# results from bayesian model


### global values ####
#load Data High
data_HIGH ### load data from read.csv or direclty from Rstudio 

####probability of death in each HIGH partition####
####Age partition####
#subset 0
#effectif in each treatment group
n01_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE >= 65, na.rm = T)
n00_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE < 65, na.rm = T)
#response in each treatment group
r01_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE >= 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r00_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE < 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)

# sous groupe 1
#effectif in each treatment group
n11_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE >= 65, na.rm = T)
n10_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE < 65, na.rm = T)
#response in each treatment group
r11_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE >= 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r10_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE < 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)

pi11_HIGH_age <- r11_age / n11_age
pi10_HIGH_age <- r10_age / n10_age
pi01_HIGH_age <- r01_age / n01_age
pi00_HIGH_age <- r00_age / n00_age
pa_HIGH_age <- mean(data_HIGH$high.AGE >= 65)

###Sofa neuro####
n01_neuro <- sum(data_HIGH$high.rando == 2 &
                   data_HIGH$high.SOFAN >= 1, na.rm = T)
n00_neuro <- sum(data_HIGH$high.rando == 2 &
                   data_HIGH$high.SOFAN < 1, na.rm = T)
#response in each treatment group
r01_neuro <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.SOFAN >= 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r00_neuro <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.SOFAN < 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)

# sous groupe 1
#effectif in each treatment group
n11_neuro <- sum(data_HIGH$high.rando == 1 &
                   data_HIGH$high.SOFAN >= 1, na.rm = T)
n10_neuro <- sum(data_HIGH$high.rando == 1 &
                   data_HIGH$high.SOFAN < 1, na.rm = T)
#response in each treatment group
r11_neuro <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.SOFAN >= 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r10_neuro <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.SOFAN < 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)

pi11_HIGH_sofa_neuro <- r11_neuro / n11_neuro
pi10_HIGH_sofa_neuro <- r10_neuro / n10_neuro
pi01_HIGH_sofa_neuro <- r01_neuro / n01_neuro
pi00_HIGH_sofa_neuro <- r00_neuro / n00_neuro
pa_HIGH_sofa_neuro <- mean(data_HIGH$high.SOFAN >= 1)

###PAO2/FIO2####
n01_respi <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.sofa_resp == 4, na.rm =
        T)
n00_respi <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.sofa_resp == 3, na.rm =
        T)
#response in each treatment group
r01_respi <-
  sum(
    data_HIGH$high.rando == 2 &
      data_HIGH$high.sofa_resp == 4 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )
r00_respi <-
  sum(
    data_HIGH$high.rando == 2 &
      data_HIGH$high.sofa_resp == 3 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )

# sous groupe 1
#effectif in each treatment group
n11_respi <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.sofa_resp == 4, na.rm =
        T)
n10_respi <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.sofa_resp == 3, na.rm =
        T)
#response in each treatment group
r11_respi <-
  sum(
    data_HIGH$high.rando == 1 &
      data_HIGH$high.sofa_resp == 4 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )
r10_respi <-
  sum(
    data_HIGH$high.rando == 1 &
      data_HIGH$high.sofa_resp == 3 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )

pi11_HIGH_sofa_resp <- r11_respi / n11_respi
pi10_HIGH_sofa_resp <- r10_respi / n10_respi
pi01_HIGH_sofa_resp <- r01_respi / n01_respi
pi00_HIGH_sofa_resp <- r00_respi / n00_respi
pa_HIGH_sofa_resp <- mean(data_HIGH$high.sofa_resp == 4)

### date of interim and terminal analysis
IA1 <- as.Date("2016-12-07", format = "%Y-%m-%d")
IA2 <- as.Date("2017-04-09", format = "%Y-%m-%d")
IA3 <- as.Date("2017-09-01", format = "%Y-%m-%d")
TA <- as.Date("2017-12-31", format = "%Y-%m-%d")

### High Binary bayesian model###
####model AGE partition####
binary_bayesian_trial_HIGH_AGE <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           eta_1,
           gamma_1,
           epsilon_1) {
    mat.sub1 <- create_matrix_HIGH()
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      Data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.AGE >= 65, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      age_65_trt <- filter(Data_HIGH_IA1, high.rando == 1 &
                             high.AGE >= 65)
      age_65_moins_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.AGE < 65)
      age_65 <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.AGE >= 65)
      age_65_moins <- filter(Data_HIGH_IA1, high.rando == 2 &
                               high.AGE < 65)
      r11_IA1 <- age_65_trt$high.dc28
      r10_IA1 <- age_65_moins_trt$high.dc28
      r01_IA1 <- age_65$high.dc28
      r00_IA1 <- age_65_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P1_0",
            "P1_1",
            "P5",
            "P6",
            "biais_0",
            "biais_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA1, Data_IA1)

      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        age_65_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.AGE < 65)
        age_65_IA2 <- filter(data_HIGH_IA2, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.AGE < 65)
        r11_IA2 <- age_65_trt_IA2$high.dc28
        r10_IA2 <- age_65_moins_trt_IA2$high.dc28
        r01_IA2 <- age_65_IA2$high.dc28
        r00_IA2 <- age_65_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
     
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_age <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                             high.Datrand <= IA2),
                                    high.AGE < 65)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 2 & high.AGE < 65)
        r10_IA2 <- age_65_moins_trt_IA2$high.dc28
        r00_IA2 <- age_65_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
     
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_age <-
          filter(filter(data_HIGH, high.Datrand > IA1 &
                          high.Datrand <= IA2),
                 high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_age)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA2_age$high.rando == 1, 1, 0)
        
        age_65_trt_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 1 &
                   high.AGE >= 65)
        
        age_65_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA2 <- age_65_trt$high.dc28
        r01_IA2 <- age_65$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
     
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        age_65_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.AGE < 65)
        age_65_IA3 <- filter(data_HIGH_IA3, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 & high.AGE < 65)
        r11_IA3 <- age_65_trt_IA3$high.dc28
        r10_IA3 <- age_65_moins_trt_IA3$high.dc28
        r01_IA3 <- age_65_IA3$high.dc28
        r00_IA3 <- age_65_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 1" |
        Check_interim_decision2[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA3_age <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                             high.Datrand <= IA3),
                                    high.AGE < 65)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 2 & high.AGE < 65)
        r10_IA3 <- age_65_moins_trt_IA3$high.dc28
        r00_IA3 <- age_65_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)

        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 0" |
        Check_interim_decision2[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA3_age <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                             high.Datrand <= IA3),
                                    high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_age$high.rando == 1, 1, 0)
        age_65_trt_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 1 &
                   high.AGE >= 65)
        age_65_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA3 <- age_65_trt_IA3$high.dc28
        r01_IA3 <- age_65_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        

        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        age_65_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.AGE < 65)
        age_65_IA4 <- filter(data_HIGH_IA4, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 & high.AGE < 65)
        r11_IA4 <- age_65_trt_IA4$high.dc28
        r10_IA4 <- age_65_moins_trt_IA4$high.dc28
        r01_IA4 <- age_65_IA4$high.dc28
        r00_IA4 <- age_65_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
 
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 1" |
        Check_interim_decision3[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA4_age <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                             high.Datrand <= TA),
                                    high.AGE < 65)
        subset <- c(rep(0, dim(data_HIGH_IA4_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 2 & high.AGE < 65)
        r10_IA4 <- age_65_moins_trt_IA4$high.dc28
        r00_IA4 <- age_65_moins_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
    
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 0" |
        Check_interim_decision3[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA4_age <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                             high.Datrand <= TA),
                                    high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_age$high.rando == 1, 1, 0)
        age_65_trt_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 1 &
                   high.AGE >= 65)
        age_65_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA4 <- age_65_trt_IA4$high.dc28
        r01_IA4 <- age_65_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
   
        Final_result <- list(Data_IA4, results_IA4)
        prin(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision4)
   
      
      mat.sub1[i, ] <-
        cbind(
          as.numeric(Final_result[[2]][1]),
          as.numeric(Final_result[[2]][2]),
          as.numeric(Final_result[[2]][3]),
          as.numeric(Final_result[[2]][4]),
          as.numeric(Final_result[[2]][5]),
          as.numeric(Final_result[[2]][6]),
          as.numeric(Final_result[[2]][7]),
          as.numeric(Final_result[[2]][8]),
          as.numeric(Final_result[[2]][9]),
          as.numeric(Check_interim_decision1[2]),
          as.numeric(Check_interim_decision1[3]),
          as.numeric(Check_interim_decision1[4]),
          as.numeric(Check_interim_decision1[5]),
          as.numeric(Check_interim_decision1[6]),
          as.numeric(Check_interim_decision1[7]),
          as.numeric(Check_interim_decision1[8]),
          as.numeric(Check_interim_decision2[2]),
          as.numeric(Check_interim_decision2[3]),
          as.numeric(Check_interim_decision2[4]),
          as.numeric(Check_interim_decision2[5]),
          as.numeric(Check_interim_decision2[6]),
          as.numeric(Check_interim_decision2[7]),
          as.numeric(Check_interim_decision2[8]),
          as.numeric(Check_interim_decision3[2]),
          as.numeric(Check_interim_decision3[3]),
          as.numeric(Check_interim_decision3[4]),
          as.numeric(Check_interim_decision3[5]),
          as.numeric(Check_interim_decision3[6]),
          as.numeric(Check_interim_decision3[7]),
          as.numeric(Check_interim_decision3[8]),
          as.numeric(Check_interim_decision4[2]),
          as.numeric(Check_interim_decision4[3]),
          as.numeric(Check_interim_decision4[4]),
          as.numeric(Check_interim_decision4[5]),
          as.numeric(Check_interim_decision4[6]),
          as.numeric(Check_interim_decision4[7]),
          as.numeric(Check_interim_decision4[8]),
          as.numeric(First_result[[2]][14]),
          as.numeric(First_result[[2]][15]),
          as.numeric(Second_result[[2]][14]),
          as.numeric(Second_result[[2]][15]),
          as.numeric(Third_result[[2]][14]),
          as.numeric(Third_result[[2]][15]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_decision[3]),
          as.numeric(Final_decision[4]),
          as.numeric(Final_decision[5]),
          as.numeric(Final_decision[6]),
          as.numeric(Final_decision[7]),
          as.numeric(Final_decision[8]),
          as.numeric(Final_decision[9]),
          as.numeric(Final_decision[10])
        )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }

#### model sofa neuro partition#####
binary_bayesian_trial_HIGH_SOFA_NEURO <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           eta_1,
           gamma_1,
           epsilon_1) {
    mat.sub1 <- create_matrix_HIGH()
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      Data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.SOFAN >= 1, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      Neuro_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.SOFAN >= 1)
      Neuro_moins_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.SOFAN < 1)
      Neuro <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.SOFAN >= 1)
      Neuro_moins <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.SOFAN < 1)
      r11_IA1 <- Neuro_trt$high.dc28
      r10_IA1 <- Neuro_moins_trt$high.dc28
      r01_IA1 <- Neuro$high.dc28
      r00_IA1 <- Neuro_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P0_0",
            "P1_0",
            "P0_1",
            "P1_1",
            "P5_0",
            "P5_1",
            "interaction_0",
            "interaction_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA1, Data_IA1)
      
 
      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        Neuro_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.SOFAN >= 1)
        Neuro_moins_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.SOFAN < 1)
        r11_IA2 <- Neuro_trt_IA2$high.dc28
        r10_IA2 <- Neuro_moins_trt_IA2$high.dc28
        r01_IA2 <- Neuro_IA2$high.dc28
        r00_IA2 <- Neuro_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
       
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_neuro <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_neuro$high.rando == 1, 1, 0)
        Neuro_moins_trt_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_moins_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r10_IA2 <- Neuro_moins_trt_IA2$high.dc28
        r00_IA2 <- Neuro_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
       
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_neuro <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_neuro)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA2_neuro$high.rando == 1, 1, 0)
        
        Neuro_trt_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        
        Neuro_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA2 <- Neuro_trt_IA2$high.dc28
        r01_IA2 <- Neuro_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
     
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        Neuro_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA3 <- filter(data_HIGH_IA3, high.rando == 2 &
                              high.SOFAN >= 1)
        Neuro_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 & high.SOFAN < 1)
        r11_IA3 <- Neuro_trt_IA3$high.dc28
        r10_IA3 <- Neuro_moins_trt_IA3$high.dc28
        r01_IA3 <- Neuro_IA3$high.dc28
        r00_IA3 <- Neuro_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
      
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 1" |
        Check_interim_decision2[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA3_neuro <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_neuro$high.rando == 1, 1, 0)
        Neuro_moins_trt_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_moins_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r10_IA3 <- Neuro_moins_trt_IA3$high.dc28
        r00_IA3 <- Neuro_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
        #glm model for gail and simon test
     
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 0" |
        Check_interim_decision2[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA3_neuro <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_neuro$high.rando == 1, 1, 0)
        Neuro_trt_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        Neuro_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA3 <- Neuro_trt_IA3$high.dc28
        r01_IA3 <- Neuro_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
        #glm model for gail and simon test
    
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        Neuro_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA4 <- filter(data_HIGH_IA4, high.rando == 2 &
                              high.SOFAN >= 1)
        Neuro_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 & high.SOFAN < 1)
        r11_IA4 <- Neuro_trt_IA4$high.dc28
        r10_IA4 <- Neuro_moins_trt_IA4$high.dc28
        r01_IA4 <- Neuro_IA4$high.dc28
        r00_IA4 <- Neuro_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        

      
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 1" |
        Check_interim_decision3[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA4_neuro <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_neuro$high.rando == 1, 1, 0)
        Neuro_trt_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r11_IA4 <- Neuro_trt_IA4$high.dc28
        r01_IA4 <- Neuro_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
      
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 0" |
        Check_interim_decision3[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA4_neuro <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_neuro$high.rando == 1, 1, 0)
        Neuro_trt_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        Neuro_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA4 <- Neuro_trt_IA4$high.dc28
        r01_IA4 <- Neuro_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
        #glm model for gail and simon test
      
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision4)
     
      
      mat.sub1[i, ] <- cbind(
        as.numeric(Final_result[[2]][1]),
        as.numeric(Final_result[[2]][2]),
        as.numeric(Final_result[[2]][3]),
        as.numeric(Final_result[[2]][4]),
        as.numeric(Final_result[[2]][5]),
        as.numeric(Final_result[[2]][6]),
        as.numeric(Final_result[[2]][7]),
        as.numeric(Final_result[[2]][8]),
        as.numeric(Final_result[[2]][9]),
        as.numeric(Check_interim_decision1[2]),
        as.numeric(Check_interim_decision1[3]),
        as.numeric(Check_interim_decision1[4]),
        as.numeric(Check_interim_decision1[5]),
        as.numeric(Check_interim_decision1[6]),
        as.numeric(Check_interim_decision1[7]),
        as.numeric(Check_interim_decision1[8]),
        as.numeric(Check_interim_decision2[2]),
        as.numeric(Check_interim_decision2[3]),
        as.numeric(Check_interim_decision2[4]),
        as.numeric(Check_interim_decision2[5]),
        as.numeric(Check_interim_decision2[6]),
        as.numeric(Check_interim_decision2[7]),
        as.numeric(Check_interim_decision2[8]),
        as.numeric(Check_interim_decision3[2]),
        as.numeric(Check_interim_decision3[3]),
        as.numeric(Check_interim_decision3[4]),
        as.numeric(Check_interim_decision3[5]),
        as.numeric(Check_interim_decision3[6]),
        as.numeric(Check_interim_decision3[7]),
        as.numeric(Check_interim_decision3[8]),
        as.numeric(Check_interim_decision4[2]),
        as.numeric(Check_interim_decision4[3]),
        as.numeric(Check_interim_decision4[4]),
        as.numeric(Check_interim_decision4[5]),
        as.numeric(Check_interim_decision4[6]),
        as.numeric(Check_interim_decision4[7]),
        as.numeric(Check_interim_decision4[8]),
        as.numeric(First_result[[2]][14]),
        as.numeric(First_result[[2]][15]),
        as.numeric(Second_result[[2]][14]),
        as.numeric(Second_result[[2]][15]),
        as.numeric(Third_result[[2]][14]),
        as.numeric(Third_result[[2]][15]),
        as.numeric(Final_result[[2]][14]),
        as.numeric(Final_result[[2]][15]),
        as.numeric(Final_decision[3]),
        as.numeric(Final_decision[4]),
        as.numeric(Final_decision[5]),
        as.numeric(Final_decision[6]),
        as.numeric(Final_decision[7]),
        as.numeric(Final_decision[8]),
        as.numeric(Final_decision[9]),
        as.numeric(Final_decision[10])
      )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }



#####pAO2/FIO2 partition #####
binary_bayesian_trial_HIGH_pao2 <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           eta_1,
           gamma_1,
           epsilon_1) {
    mat.sub1 <- create_matrix_HIGH()
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.sofa_resp == 4, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      PAO2_trt <-
        filter(data_HIGH_IA1, high.rando == 1 & high.sofa_resp == 4)
      PAO2_moins_trt <-
        filter(data_HIGH_IA1, high.rando == 1 & high.sofa_resp == 3)
      PAO2 <-
        filter(data_HIGH_IA1, high.rando == 2 & high.sofa_resp == 4)
      PAO2_moins <-
        filter(data_HIGH_IA1, high.rando == 2 & high.sofa_resp == 3)
      r11_IA1 <- PAO2_trt$high.dc28
      r10_IA1 <- PAO2_moins_trt$high.dc28
      r01_IA1 <- PAO2$high.dc28
      r00_IA1 <- PAO2_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P0_0",
            "P1_0",
            "P0_1",
            "P1_1",
            "P5_0",
            "P5_1",
            "interaction_0",
            "interaction_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA1, Data_IA1)
      
 
      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        PAO2_trt <-
          filter(data_HIGH_IA2, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt <-
          filter(data_HIGH_IA2, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2 <-
          filter(data_HIGH_IA2, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins <-
          filter(data_HIGH_IA2, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA2 <- PAO2_trt_IA2$high.dc28
        r10_IA2 <- PAO2_moins_trt_IA2$high.dc28
        r01_IA2 <- PAO2_IA2$high.dc28
        r00_IA2 <- PAO2_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
       
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_respi <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.sofa_resp == 3)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA2 <- PAO2_moins_trt_IA2$high.dc28
        r00_IA2 <- PAO2_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
       
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_respi <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.sofa_resp == 4)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_respi)[1]))
        group <- NULL
        group <- group <-
          ifelse(data_HIGH_IA2_respi$high.rando == 1, 1, 0)
        
        PAO2_trt_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        
        PAO2_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA2 <- PAO2_trt$high.dc28
        r01_IA2 <- PAO2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA2, Data_IA2)
        
      
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        PAO2_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA3 <- PAO2_trt_IA3$high.dc28
        r10_IA3 <- PAO2_moins_trt_IA3$high.dc28
        r01_IA3 <- PAO2_IA3$high.dc28
        r00_IA3 <- PAO2_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
     
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 1" |
        Check_interim_decision2[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA3_respi <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.sofa_resp == 3)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA3 <- PAO2_moins_trt_IA3$high.dc28
        r00_IA3 <- PAO2_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
      
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 0" |
        Check_interim_decision2[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA3_respi <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.sofa_resp == 4)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_respi$high.rando == 1, 1, 0)
        PAO2_trt_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA3 <- PAO2_trt_IA3$high.dc28
        r01_IA3 <- PAO2_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA3, Data_IA3)
        
      
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        PAO2_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA4 <- PAO2_trt_IA4$high.dc28
        r10_IA4 <- PAO2_moins_trt_IA4$high.dc28
        r01_IA4 <- PAO2_IA4$high.dc28
        r00_IA4 <- PAO2_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
      
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 1" |
        Check_interim_decision3[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA4_respi <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.sofa_resp == 3)
        subset <- c(rep(0, dim(data_HIGH_IA4_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA4 <- PAO2_moins_trt_IA4$high.dc28
        r00_IA4 <- PAO2_moins_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
     
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 0" |
        Check_interim_decision3[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA4_respi <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.sofa_resp == 4)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_respi$high.rando == 1, 1, 0)
        PAO2_trt_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA4 <- PAO2_trt_IA4$high.dc28
        r01_IA4 <- PAO2_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P0_0",
              "P1_0",
              "P0_1",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_HIGH_DATA(bayesian_jags_fonction_IA4, Data_IA4)
        
     
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision4)
    
      
      mat.sub1[i, ] <-
        cbind(
          as.numeric(Final_result[[2]][1]),
          as.numeric(Final_result[[2]][2]),
          as.numeric(Final_result[[2]][3]),
          as.numeric(Final_result[[2]][4]),
          as.numeric(Final_result[[2]][5]),
          as.numeric(Final_result[[2]][6]),
          as.numeric(Final_result[[2]][7]),
          as.numeric(Final_result[[2]][8]),
          as.numeric(Final_result[[2]][9]),
          as.numeric(Check_interim_decision1[2]),
          as.numeric(Check_interim_decision1[3]),
          as.numeric(Check_interim_decision1[4]),
          as.numeric(Check_interim_decision1[5]),
          as.numeric(Check_interim_decision1[6]),
          as.numeric(Check_interim_decision1[7]),
          as.numeric(Check_interim_decision1[8]),
          as.numeric(Check_interim_decision2[2]),
          as.numeric(Check_interim_decision2[3]),
          as.numeric(Check_interim_decision2[4]),
          as.numeric(Check_interim_decision2[5]),
          as.numeric(Check_interim_decision2[6]),
          as.numeric(Check_interim_decision2[7]),
          as.numeric(Check_interim_decision2[8]),
          as.numeric(Check_interim_decision3[2]),
          as.numeric(Check_interim_decision3[3]),
          as.numeric(Check_interim_decision3[4]),
          as.numeric(Check_interim_decision3[5]),
          as.numeric(Check_interim_decision3[6]),
          as.numeric(Check_interim_decision3[7]),
          as.numeric(Check_interim_decision3[8]),
          as.numeric(Check_interim_decision4[2]),
          as.numeric(Check_interim_decision4[3]),
          as.numeric(Check_interim_decision4[4]),
          as.numeric(Check_interim_decision4[5]),
          as.numeric(Check_interim_decision4[6]),
          as.numeric(Check_interim_decision4[7]),
          as.numeric(Check_interim_decision4[8]),
          as.numeric(First_result[[2]][14]),
          as.numeric(First_result[[2]][15]),
          as.numeric(Second_result[[2]][14]),
          as.numeric(Second_result[[2]][15]),
          as.numeric(Third_result[[2]][14]),
          as.numeric(Third_result[[2]][15]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_decision[3]),
          as.numeric(Final_decision[4]),
          as.numeric(Final_decision[5]),
          as.numeric(Final_decision[6]),
          as.numeric(Final_decision[7]),
          as.numeric(Final_decision[8]),
          as.numeric(Final_decision[9]),
          as.numeric(Final_decision[10])
        )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }




##### results from millen method ####
Nsimu = 1

set.seed(18121995)

neuro_millen <- binary_bayesian_trial_HIGH_SOFA_NEURO(
  pi11_HIGH_sofa_neuro,
  pi10_HIGH_sofa_neuro,
  pi01_HIGH_sofa_neuro,
  pi00_HIGH_sofa_neuro,
  pa_HIGH_sofa_neuro,
  binary_bayesian_model_2_subset_Millen,
  stopping_rules_first_analyse_variation_2_subset_millen,
  stopping_rules_other_analyse_variation_2_subset_millen,
  0.9,
  1.25,
  0.90,
  0.90
)
age_millen <- binary_bayesian_trial_HIGH_AGE(
  pi11_HIGH_age,
  pi10_HIGH_age,
  pi01_HIGH_age,
  pi00_HIGH_age,
  pa_HIGH_age,
  binary_bayesian_model_2_subset_Millen,
  stopping_rules_first_analyse_variation_2_subset_millen,
  stopping_rules_other_analyse_variation_2_subset_millen,
  0.9,1.25,0.90,0.85
)

paO2_millen <- binary_bayesian_trial_HIGH_pao2(
  pi11_HIGH_sofa_resp,
  pi10_HIGH_sofa_resp,
  pi01_HIGH_sofa_resp,
  pi00_HIGH_sofa_resp,
  pa_HIGH_sofa_resp,
  binary_bayesian_model_2_subset_Millen,
  stopping_rules_first_analyse_variation_2_subset_millen,
  stopping_rules_other_analyse_variation_2_subset_millen,
  lambda_1,
  eta_1,
  gamma_1,
  epsilon_1
)


#### Method with Gail and Simon test : load Simulations_method_with GS_interaction for all functions ####
#####age#####
binary_bayesian_trial_HIGH_AGE_GS <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           gamma_1,
           epsilon_1,
           C1,
           C2) {
    mat.sub1 <- create_matrix_2.0(Nsimu)
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      Data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.AGE >= 65, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      age_65_trt <- filter(Data_HIGH_IA1, high.rando == 1 &
                             high.AGE >= 65)
      age_65_moins_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.AGE < 65)
      age_65 <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.AGE >= 65)
      age_65_moins <- filter(Data_HIGH_IA1, high.rando == 2 &
                               high.AGE < 65)
      r11_IA1 <- age_65_trt$high.dc28
      r10_IA1 <- age_65_moins_trt$high.dc28
      r01_IA1 <- age_65$high.dc28
      r00_IA1 <- age_65_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values_2.0(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P1_0",
            "P1_1",
            "P5",
            "P6",
            "biais_0",
            "biais_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_2.0(bayesian_jags_fonction_IA1, Data_IA1)
      
      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        age_65_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.AGE < 65)
        age_65_IA2 <- filter(data_HIGH_IA2, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.AGE < 65)
        r11_IA2 <- age_65_trt_IA2$high.dc28
        r10_IA2 <- age_65_moins_trt_IA2$high.dc28
        r01_IA2 <- age_65_IA2$high.dc28
        r00_IA2 <- age_65_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_age <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                             high.Datrand <= IA2),
                                    high.AGE < 65)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 2 & high.AGE < 65)
        r10_IA2 <- age_65_moins_trt_IA2$high.dc28
        r00_IA2 <- age_65_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_age <-
          filter(filter(data_HIGH, high.Datrand > IA1 &
                          high.Datrand <= IA2),
                 high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_age)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA2_age$high.rando == 1, 1, 0)
        
        age_65_trt_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 1 &
                   high.AGE >= 65)
        
        age_65_IA2 <-
          filter(data_HIGH_IA2_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA2 <- age_65_trt$high.dc28
        r01_IA2 <- age_65$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision1[9]
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        age_65_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.AGE < 65)
        age_65_IA3 <- filter(data_HIGH_IA3, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 & high.AGE < 65)
        r11_IA3 <- age_65_trt_IA3$high.dc28
        r10_IA3 <- age_65_moins_trt_IA3$high.dc28
        r01_IA3 <- age_65_IA3$high.dc28
        r00_IA3 <- age_65_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 1" |
        Check_interim_decision2[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA3_age <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                             high.Datrand <= IA3),
                                    high.AGE < 65)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 2 & high.AGE < 65)
        r10_IA3 <- age_65_moins_trt_IA3$high.dc28
        r00_IA3 <- age_65_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "futility in the subset 0" |
        Check_interim_decision2[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA3_age <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                             high.Datrand <= IA3),
                                    high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_age$high.rando == 1, 1, 0)
        age_65_trt_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 1 &
                   high.AGE >= 65)
        age_65_IA3 <-
          filter(data_HIGH_IA3_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA3 <- age_65_trt_IA3$high.dc28
        r01_IA3 <- age_65_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        #glm model for gail and simon test
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision2[9]
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.AGE >= 65, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        age_65_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.AGE >= 65)
        age_65_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.AGE < 65)
        age_65_IA4 <- filter(data_HIGH_IA4, high.rando == 2 &
                               high.AGE >= 65)
        age_65_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 & high.AGE < 65)
        r11_IA4 <- age_65_trt_IA4$high.dc28
        r10_IA4 <- age_65_moins_trt_IA4$high.dc28
        r01_IA4 <- age_65_IA4$high.dc28
        r00_IA4 <- age_65_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 1" |
        Check_interim_decision3[1] == "efficacy in the subset 0"
      )) {
        data_HIGH_IA4_age <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                             high.Datrand <= TA),
                                    high.AGE < 65)
        subset <- c(rep(0, dim(data_HIGH_IA4_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_age$high.rando == 1, 1, 0)
        age_65_moins_trt_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 1 & high.AGE < 65)
        age_65_moins_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 2 & high.AGE < 65)
        r10_IA4 <- age_65_moins_trt_IA4$high.dc28
        r00_IA4 <- age_65_moins_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "futility in the subset 0" |
        Check_interim_decision3[1] == "efficacy in the subset 1"
      )) {
        data_HIGH_IA4_age <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                             high.Datrand <= TA),
                                    high.AGE >= 65)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_age)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_age$high.rando == 1, 1, 0)
        age_65_trt_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 1 &
                   high.AGE >= 65)
        age_65_IA4 <-
          filter(data_HIGH_IA4_age, high.rando == 2 &
                   high.AGE >= 65)
        r11_IA4 <- age_65_trt_IA4$high.dc28
        r01_IA4 <- age_65_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        prin(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision3[9]
        )
      print(Check_interim_decision4)
      Final_decision <-
        stopping_rules_Terminal_analyse(Final_result[[2]][4], Final_result[[2]][7])
      print(Final_decision)
      
      mat.sub1 <-
        cbind(
          gamma,
          C1,
          C2,
          pa,
          as.numeric(Final_result[[2]][1]),
          as.numeric(Final_result[[2]][2]),
          as.numeric(Final_result[[2]][3]),
          as.numeric(Final_result[[2]][4]),
          as.numeric(Final_result[[2]][5]),
          as.numeric(Final_result[[2]][6]),
          as.numeric(Final_result[[2]][7]),
          as.numeric(Final_result[[2]][8]),
          as.numeric(Final_result[[2]][9]),
          as.numeric(Check_interim_decision1[2]),
          as.numeric(Check_interim_decision1[3]),
          as.numeric(Check_interim_decision1[4]),
          as.numeric(Check_interim_decision1[5]),
          as.numeric(Check_interim_decision1[6]),
          as.numeric(Check_interim_decision1[7]),
          as.numeric(Check_interim_decision1[8]),
          as.numeric(Check_interim_decision2[2]),
          as.numeric(Check_interim_decision2[3]),
          as.numeric(Check_interim_decision2[4]),
          as.numeric(Check_interim_decision2[5]),
          as.numeric(Check_interim_decision2[6]),
          as.numeric(Check_interim_decision2[7]),
          as.numeric(Check_interim_decision2[8]),
          as.numeric(Check_interim_decision3[2]),
          as.numeric(Check_interim_decision3[3]),
          as.numeric(Check_interim_decision3[4]),
          as.numeric(Check_interim_decision3[5]),
          as.numeric(Check_interim_decision3[6]),
          as.numeric(Check_interim_decision3[7]),
          as.numeric(Check_interim_decision3[8]),
          as.numeric(Check_interim_decision4[2]),
          as.numeric(Check_interim_decision4[3]),
          as.numeric(Check_interim_decision4[4]),
          as.numeric(Check_interim_decision4[5]),
          as.numeric(Check_interim_decision4[6]),
          as.numeric(Check_interim_decision4[7]),
          as.numeric(Check_interim_decision4[8]),
          as.numeric(First_result[[2]][18]),
          as.numeric(First_result[[2]][19]),
          as.numeric(Second_result[[2]][18]),
          as.numeric(Second_result[[2]][19]),
          as.numeric(Third_result[[2]][18]),
          as.numeric(Third_result[[2]][19]),
          as.numeric(Final_result[[2]][18]),
          as.numeric(Final_result[[2]][19]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_result[[2]][16]),
          as.numeric(Final_result[[2]][17])
        )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }


age_GS <- binary_bayesian_trial_HIGH_AGE_GS(
  pi11_HIGH_age,
  pi10_HIGH_age,
  pi01_HIGH_age,
  pi00_HIGH_age,
  pa_HIGH_age,
  binary_bayesian_model_2_subset_GS,
  stopping_rules_first_analyse_variation_2_subset_GS,
  stopping_rules_other_analyse_variation_2_subset_GS,
  0.9,
  0.85,
  0.10,
  0.465,
  0.215
)

#####sofa neuro####
binary_bayesian_trial_HIGH_sofa_neuro_GS <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           gamma_1,
           epsilon_1,
           C1,
           C2) {
    mat.sub1 <- create_matrix_2.0(Nsimu)
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      Data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.SOFAN >= 1, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      Neuro_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.SOFAN >= 1)
      Neuro_moins_trt <-
        filter(Data_HIGH_IA1, high.rando == 1 & high.SOFAN < 1)
      Neuro <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.SOFAN >= 1)
      Neuro_moins <-
        filter(Data_HIGH_IA1, high.rando == 2 & high.SOFAN < 1)
      r11_IA1 <- Neuro_trt$high.dc28
      r10_IA1 <- Neuro_moins_trt$high.dc28
      r01_IA1 <- Neuro$high.dc28
      r00_IA1 <- Neuro_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values_2.0(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P1_0",
            "P1_1",
            "P5",
            "P6",
            "biais_0",
            "biais_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_2.0(bayesian_jags_fonction_IA1, Data_IA1)
      
      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        Neuro_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA2 <-
          filter(data_HIGH_IA2, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.SOFAN >= 1)
        Neuro_moins_IA2 <-
          filter(data_HIGH_IA2, high.rando == 2 & high.SOFAN < 1)
        r11_IA2 <- Neuro_trt_IA2$high.dc28
        r10_IA2 <- Neuro_moins_trt_IA2$high.dc28
        r01_IA2 <- Neuro_IA2$high.dc28
        r00_IA2 <- Neuro_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_neuro <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_neuro$high.rando == 1, 1, 0)
        Neuro_moins_trt_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_moins_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r10_IA2 <- Neuro_moins_trt_IA2$high.dc28
        r00_IA2 <- Neuro_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_neuro <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_neuro)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA2_neuro$high.rando == 1, 1, 0)
        
        Neuro_trt_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        
        Neuro_IA2 <-
          filter(data_HIGH_IA2_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA2 <- Neuro_trt_IA2$high.dc28
        r01_IA2 <- Neuro_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision1[9]
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        Neuro_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 & high.SOFAN >= 1)
        Neuro_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 & high.SOFAN < 1)
        r11_IA3 <- Neuro_trt_IA3$high.dc28
        r10_IA3 <- Neuro_moins_trt_IA3$high.dc28
        r01_IA3 <- Neuro_IA3$high.dc28
        r00_IA3 <- Neuro_moins_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if (Check_interim_decision2[1] == "futility in the subset 1" |
                 Check_interim_decision2[1] == "efficacy in the subset 0") {
        data_HIGH_IA3_neuro <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_neuro$high.rando == 1, 1, 0)
        Neuro_moins_trt_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_moins_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r10_IA3 <- Neuro_moins_trt_IA3$high.dc28
        r00_IA3 <- Neuro_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      } else if (Check_interim_decision2[1] == "futility in the subset 0" |
                 Check_interim_decision2[1] == "efficacy in the subset 1") {
        data_HIGH_IA3_neuro <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_neuro)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA3_neuro$high.rando == 1, 1, 0)
        
        Neuro_trt_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        
        Neuro_IA3 <-
          filter(data_HIGH_IA3_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA3 <- Neuro_trt_IA3$high.dc28
        r01_IA3 <- Neuro_IA3$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        Data_IA3 <- data.frame(group, subset, y)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision2[9]
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.SOFAN >= 1, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        Neuro_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.SOFAN >= 1)
        Neuro_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 & high.SOFAN < 1)
        Neuro_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 & high.SOFAN >= 1)
        Neuro_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 & high.SOFAN < 1)
        r11_IA4 <- Neuro_trt_IA4$high.dc28
        r10_IA4 <- Neuro_moins_trt_IA4$high.dc28
        r01_IA4 <- Neuro_IA4$high.dc28
        r00_IA4 <- Neuro_moins_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
        
        
      } else if (Check_interim_decision3[1] == "futility in the subset 1" |
                 Check_interim_decision3[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_neuro <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.SOFAN < 1)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA4_neuro)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_neuro$high.rando == 1, 1, 0)
        Neuro_moins_trt_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 1 &
                   high.SOFAN < 1)
        Neuro_moins_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 2 &
                   high.SOFAN < 1)
        r10_IA4 <- Neuro_moins_trt_IA4$high.dc28
        r00_IA4 <- Neuro_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      } else if (Check_interim_decision3[1] == "futility in the subset 0" |
                 Check_interim_decision3[1] == "efficacy in the subset 1") {
        data_HIGH_IA4_neuro <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.SOFAN >= 1)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA4_neuro)[1]))
        group <- NULL
        group <-
          group <- ifelse(data_HIGH_IA4_neuro$high.rando == 1, 1, 0)
        
        Neuro_trt_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 1 &
                   high.SOFAN >= 1)
        
        Neuro_IA4 <-
          filter(data_HIGH_IA4_neuro, high.rando == 2 &
                   high.SOFAN >= 1)
        r11_IA4 <- Neuro_trt_IA4$high.dc28
        r01_IA4 <- Neuro_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision3[9]
        )
      print(Check_interim_decision4)
      Final_decision <-
        stopping_rules_Terminal_analyse(Final_result[[2]][4], Final_result[[2]][7])
      print(Final_decision)
      
      mat.sub1 <-
        cbind(
          gamma,
          C1,
          C2,
          pa,
          as.numeric(Final_result[[2]][1]),
          as.numeric(Final_result[[2]][2]),
          as.numeric(Final_result[[2]][3]),
          as.numeric(Final_result[[2]][4]),
          as.numeric(Final_result[[2]][5]),
          as.numeric(Final_result[[2]][6]),
          as.numeric(Final_result[[2]][7]),
          as.numeric(Final_result[[2]][8]),
          as.numeric(Final_result[[2]][9]),
          as.numeric(Check_interim_decision1[2]),
          as.numeric(Check_interim_decision1[3]),
          as.numeric(Check_interim_decision1[4]),
          as.numeric(Check_interim_decision1[5]),
          as.numeric(Check_interim_decision1[6]),
          as.numeric(Check_interim_decision1[7]),
          as.numeric(Check_interim_decision1[8]),
          as.numeric(Check_interim_decision2[2]),
          as.numeric(Check_interim_decision2[3]),
          as.numeric(Check_interim_decision2[4]),
          as.numeric(Check_interim_decision2[5]),
          as.numeric(Check_interim_decision2[6]),
          as.numeric(Check_interim_decision2[7]),
          as.numeric(Check_interim_decision2[8]),
          as.numeric(Check_interim_decision3[2]),
          as.numeric(Check_interim_decision3[3]),
          as.numeric(Check_interim_decision3[4]),
          as.numeric(Check_interim_decision3[5]),
          as.numeric(Check_interim_decision3[6]),
          as.numeric(Check_interim_decision3[7]),
          as.numeric(Check_interim_decision3[8]),
          as.numeric(Check_interim_decision4[2]),
          as.numeric(Check_interim_decision4[3]),
          as.numeric(Check_interim_decision4[4]),
          as.numeric(Check_interim_decision4[5]),
          as.numeric(Check_interim_decision4[6]),
          as.numeric(Check_interim_decision4[7]),
          as.numeric(Check_interim_decision4[8]),
          as.numeric(First_result[[2]][18]),
          as.numeric(First_result[[2]][19]),
          as.numeric(Second_result[[2]][18]),
          as.numeric(Second_result[[2]][19]),
          as.numeric(Third_result[[2]][18]),
          as.numeric(Third_result[[2]][19]),
          as.numeric(Final_result[[2]][18]),
          as.numeric(Final_result[[2]][19]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_result[[2]][16]),
          as.numeric(Final_result[[2]][17])
        )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }
neuro_GS <- binary_bayesian_trial_HIGH_sofa_neuro_GS(
  pi11_HIGH_sofa_neuro,
  pi10_HIGH_sofa_neuro,
  pi01_HIGH_sofa_neuro,
  pi00_HIGH_sofa_neuro,
  pa_HIGH_sofa_neuro,
  binary_bayesian_model_2_subset_GS,
  stopping_rules_first_analyse_variation_2_subset_GS,
  stopping_rules_other_analyse_variation_2_subset_GS,
  0.90,
  0.90,
  0.05,
  0.655,
  0.220
)



#### POA2/FIO2#####
binary_bayesian_trial_HIGH_sofa_resp_GS <-
  function(pi11,
           pi10,
           pi01,
           pi00,
           pa,
           bayesian.model,
           stopping_rules.IA,
           stopping_rules.TA,
           lambda_1,
           gamma_1,
           epsilon_1,
           C1,
           C2) {
    mat.sub1 <- create_matrix_2.0(Nsimu)
    
    for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      ##### HIGH data management####
      data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
      subset <- NULL
      subset <- ifelse(Data_HIGH_IA1$high.sofa_resp == 4, 1, 0)
      group <- NULL
      group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
      PAO2_trt <-
        filter(data_HIGH_IA1, high.rando == 1 & high.sofa_resp == 4)
      PAO2_moins_trt <-
        filter(data_HIGH_IA1, high.rando == 1 & high.sofa_resp == 3)
      PAO2 <-
        filter(data_HIGH_IA1, high.rando == 2 & high.sofa_resp == 4)
      PAO2_moins <-
        filter(data_HIGH_IA1, high.rando == 2 & high.sofa_resp == 3)
      r11_IA1 <- PAO2_trt$high.dc28
      r10_IA1 <- PAO2_moins_trt$high.dc28
      r01_IA1 <- PAO2$high.dc28
      r00_IA1 <- PAO2_moins$high.dc28
      ##### Subset biomarkers
      y <- NULL                #### Response
      y[group == 1 & subset == 1] <- r11_IA1
      y[group == 1 & subset == 0] <- r10_IA1
      y[group == 0 & subset == 1] <- r01_IA1
      y[group == 0 & subset == 0] <- r00_IA1
      Data_IA1 <- data.frame(group, subset, y)
      
      inits_IA1 <-
        create_list_for_inits_values_2.0(Data_IA1, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = binary_bayesian_model_2_subset_Millen_2.0,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim",
            "p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0",
            "RR",
            "P1_0",
            "P1_1",
            "P5",
            "P6",
            "biais_0",
            "biais_1"
          ),
          jags.seed = 18121995
        )
      
      results_IA1 <-
        results_bayesian_model_2.0(bayesian_jags_fonction_IA1, Data_IA1)
      
      First_result <- list(Data_IA1, results_IA1)
      print(results_IA1)
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],
          gamma_1,
          epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets") {
        data_HIGH_IA2 <- filter(data_HIGH, high.Datrand <= IA2)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA2$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA2$high.rando == 1, 1, 0)
        PAO2_trt <-
          filter(data_HIGH_IA2, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt <-
          filter(data_HIGH_IA2, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2 <-
          filter(data_HIGH_IA2, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins <-
          filter(data_HIGH_IA2, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA2 <- PAO2_trt_IA2$high.dc28
        r10_IA2 <- PAO2_moins_trt_IA2$high.dc28
        r01_IA2 <- PAO2_IA2$high.dc28
        r00_IA2 <- PAO2_moins_IA2$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA2_respi <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.sofa_resp == 3)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA2_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA2 <- PAO2_moins_trt_IA2$high.dc28
        r00_IA2 <- PAO2_moins_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA2
        y[group == 0 & subset == 0] <- r00_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA2_respi <- filter(filter(data_HIGH, high.Datrand > IA1 &
                                               high.Datrand <= IA2),
                                      high.sofa_resp == 4)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA2_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA2_respi$high.rando == 1, 1, 0)
        PAO2_trt_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_IA2 <-
          filter(data_HIGH_IA2_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA2 <- PAO2_trt_IA2$high.dc28
        r01_IA2 <- PAO2_IA2$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA2
        y[group == 0 & subset == 1] <- r01_IA2
        Data_IA2 <- data.frame(group, subset, y)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        inits_IA2 <-
          create_list_for_inits_values_2.0(Data_IA2, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA2 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA2, Data_IA2)
        
        
        Second_result <- list(Data_IA2, results_IA2)
        print(results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision1[9]
        )
      print(Check_interim_decision2)
      if ((Check_interim_decision2[1] == "continue with 2 subsets")) {
        data_HIGH_IA3 <- filter(data_HIGH, high.Datrand <= IA3)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA3$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA3$high.rando == 1, 1, 0)
        PAO2_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt_IA3 <-
          filter(data_HIGH_IA3, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins_IA3 <-
          filter(data_HIGH_IA3, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA3 <- PAO2_trt_IA3$high.dc28
        r10_IA3 <- PAO2_moins_trt_IA3$high.dc28
        r01_IA3 <- PAO2_IA3$high.dc28
        r00_IA3 <- PAO2_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        Data_IA3 <- data.frame(group, subset, y)
        
        
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA3_respi <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.sofa_resp == 3)
        subset <- NULL
        subset <- c(rep(0, dim(data_HIGH_IA3_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA3 <- PAO2_moins_trt_IA3$high.dc28
        r00_IA3 <- PAO2_moins_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA3
        y[group == 0 & subset == 0] <- r00_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA3_respi <- filter(filter(data_HIGH, high.Datrand > IA2 &
                                               high.Datrand <= IA3),
                                      high.sofa_resp == 4)
        subset <- NULL
        subset <- c(rep(1, dim(data_HIGH_IA3_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA3_respi$high.rando == 1, 1, 0)
        PAO2_trt_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_IA3 <-
          filter(data_HIGH_IA3_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA3 <- PAO2_trt_IA3$high.dc28
        r01_IA3 <- PAO2_IA3$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA3
        y[group == 0 & subset == 1] <- r01_IA3
        
        Data_IA3 <- data.frame(group, subset, y)
        inits_IA3 <-
          create_list_for_inits_values_2.0(Data_IA3, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA3 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA3, Data_IA3)
        
        
        Third_result <- list(Data_IA3, results_IA3)
        print(results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision2[9]
        )
      print(Check_interim_decision3)
      if ((Check_interim_decision3[1] == "continue with 2 subsets")) {
        data_HIGH_IA4 <- filter(data_HIGH, high.Datrand <= TA)
        subset <- NULL
        subset <- ifelse(data_HIGH_IA4$high.sofa_resp == 4, 1, 0)
        group <- NULL
        group <- ifelse(data_HIGH_IA4$high.rando == 1, 1, 0)
        PAO2_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_moins_trt_IA4 <-
          filter(data_HIGH_IA4, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 &
                   high.sofa_resp == 4)
        PAO2_moins_IA4 <-
          filter(data_HIGH_IA4, high.rando == 2 &
                   high.sofa_resp == 3)
        r11_IA4 <- PAO2_trt_IA4$high.dc28
        r10_IA4 <- PAO2_moins_trt_IA4$high.dc28
        r01_IA4 <- PAO2_IA4$high.dc28
        r00_IA4 <- PAO2_moins_IA4$high.dc28
        
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "futility in the subset 1" |
                 Check_interim_decision1[1] == "efficacy in the subset 0") {
        data_HIGH_IA4_respi <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.sofa_resp == 3)
        subset <- c(rep(0, dim(data_HIGH_IA4_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_respi$high.rando == 1, 1, 0)
        PAO2_moins_trt_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 1 &
                   high.sofa_resp == 3)
        PAO2_moins_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 2 &
                   high.sofa_resp == 3)
        r10_IA4 <- PAO2_moins_trt_IA4$high.dc28
        r00_IA4 <- PAO2_moins_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 0] <- r10_IA4
        y[group == 0 & subset == 0] <- r00_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      } else if (Check_interim_decision1[1] == "futility in the subset 0" |
                 Check_interim_decision1[1] == "efficacy in the subset 1") {
        data_HIGH_IA4_respi <- filter(filter(data_HIGH, high.Datrand > IA3 &
                                               high.Datrand <= TA),
                                      high.sofa_resp == 4)
        subset <- c(rep(0, dim(data_HIGH_IA4_respi)[1]))
        group <- NULL
        group <- ifelse(data_HIGH_IA4_respi$high.rando == 1, 1, 0)
        PAO2_trt_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 1 &
                   high.sofa_resp == 4)
        PAO2_IA4 <-
          filter(data_HIGH_IA4_respi, high.rando == 2 &
                   high.sofa_resp == 4)
        r11_IA4 <- PAO2_trt_IA4$high.dc28
        r01_IA4 <- PAO2_IA4$high.dc28
        ##### Subset biomarkers
        y <- NULL                #### Response
        y[group == 1 & subset == 1] <- r11_IA4
        y[group == 0 & subset == 1] <- r01_IA4
        Data_IA4 <- data.frame(group, subset, y)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        inits_IA4 <-
          create_list_for_inits_values_2.0(Data_IA4, pi11, pi10, pi01, pi00, pa, lambda_1, C1, C2)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim",
              "p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0",
              "RR",
              "P1_0",
              "P1_1",
              "P5",
              "P6",
              "biais_0",
              "biais_1"
            ),
            jags.seed = 18121995
          )
        
        results_IA4 <-
          results_bayesian_model_2.0(bayesian_jags_fonction_IA4, Data_IA4)
        
        
        Final_result <- list(Data_IA4, results_IA4)
        print(results_IA4)
      }
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],
          gamma_1,
          epsilon_1,
          Check_interim_decision3[9]
        )
      print(Check_interim_decision4)
      Final_decision <-
        stopping_rules_Terminal_analyse(Final_result[[2]][4], Final_result[[2]][7])
      print(Final_decision)
      
      mat.sub1 <-
        cbind(
          gamma,
          C1,
          C2,
          q1,
          pa,
          as.numeric(Final_result[[2]][1]),
          as.numeric(Final_result[[2]][2]),
          as.numeric(Final_result[[2]][3]),
          as.numeric(Final_result[[2]][4]),
          as.numeric(Final_result[[2]][5]),
          as.numeric(Final_result[[2]][6]),
          as.numeric(Final_result[[2]][7]),
          as.numeric(Final_result[[2]][8]),
          as.numeric(Final_result[[2]][9]),
          as.numeric(Check_interim_decision1[2]),
          as.numeric(Check_interim_decision1[3]),
          as.numeric(Check_interim_decision1[4]),
          as.numeric(Check_interim_decision1[5]),
          as.numeric(Check_interim_decision1[6]),
          as.numeric(Check_interim_decision1[7]),
          as.numeric(Check_interim_decision1[8]),
          as.numeric(Check_interim_decision2[2]),
          as.numeric(Check_interim_decision2[3]),
          as.numeric(Check_interim_decision2[4]),
          as.numeric(Check_interim_decision2[5]),
          as.numeric(Check_interim_decision2[6]),
          as.numeric(Check_interim_decision2[7]),
          as.numeric(Check_interim_decision2[8]),
          as.numeric(Check_interim_decision3[2]),
          as.numeric(Check_interim_decision3[3]),
          as.numeric(Check_interim_decision3[4]),
          as.numeric(Check_interim_decision3[5]),
          as.numeric(Check_interim_decision3[6]),
          as.numeric(Check_interim_decision3[7]),
          as.numeric(Check_interim_decision3[8]),
          as.numeric(Check_interim_decision4[2]),
          as.numeric(Check_interim_decision4[3]),
          as.numeric(Check_interim_decision4[4]),
          as.numeric(Check_interim_decision4[5]),
          as.numeric(Check_interim_decision4[6]),
          as.numeric(Check_interim_decision4[7]),
          as.numeric(Check_interim_decision4[8]),
          as.numeric(First_result[[2]][18]),
          as.numeric(First_result[[2]][19]),
          as.numeric(Second_result[[2]][18]),
          as.numeric(Second_result[[2]][19]),
          as.numeric(Third_result[[2]][18]),
          as.numeric(Third_result[[2]][19]),
          as.numeric(Final_result[[2]][18]),
          as.numeric(Final_result[[2]][19]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_result[[2]][16]),
          as.numeric(Final_result[[2]][17])
        )
    }
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }

resp_GS <- binary_bayesian_trial_HIGH_sofa_resp_GS(
  pi11_HIGH_sofa_resp,
  pi10_HIGH_sofa_resp,
  pi01_HIGH_sofa_resp,
  pi00_HIGH_sofa_resp,
  pa_HIGH_sofa_resp,
  binary_bayesian_model_2_subset_GS,
  stopping_rules_first_analyse_variation_2_subset_GS,
  stopping_rules_other_analyse_variation_2_subset_GS,
  lambda_1,
  0.80,
  epsilon_1,
  C1,
  C2
)



#### parition 3 sous- groupe Age  coag pas possible effectif trop différents######

 table(cut(data_HIGH$high.AGE,c(18,58,68,90),include.lowest = TRUE))
 
 n11_age <- sum(data_HIGH$high.rando == 1 &
                  data_HIGH$high.AGE <=58, na.rm = T)
 n01_age <- sum(data_HIGH$high.rando == 2 &
                  data_HIGH$high.AGE <=58, na.rm = T)
 #response in each treatment group
 r11_age <-
   sum(data_HIGH$high.rando == 1 &
         data_HIGH$high.AGE <=58 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 r01_age <-
   sum(data_HIGH$high.rando == 2 &
         data_HIGH$high.AGE <=58 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 
 # sous groupe 1
 #effectif in each treatment group
 n12_age <- sum(data_HIGH$high.rando == 1 & data_HIGH$high.AGE >58 &
                  data_HIGH$high.AGE <= 68, na.rm = T)
 n02_age <- sum(data_HIGH$high.rando == 2 & data_HIGH$high.AGE >58 &
                  data_HIGH$high.AGE <= 68, na.rm = T)
 #response in each treatment group
 r12_age <-
   sum(data_HIGH$high.rando == 1 & data_HIGH$high.AGE >58 &
         data_HIGH$high.AGE <= 68 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 r02_age <-
   sum(data_HIGH$high.rando == 2 &data_HIGH$high.AGE >58&
         data_HIGH$high.AGE <= 68 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 
 n13_age <- sum(data_HIGH$high.rando == 1 &
                  data_HIGH$high.AGE > 68, na.rm = T)
 n03_age <- sum(data_HIGH$high.rando == 2 &
                  data_HIGH$high.AGE > 68, na.rm = T)
 #response in each treatment group
 r13_age <-
   sum(data_HIGH$high.rando == 1 &
         data_HIGH$high.AGE > 68 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 r03_age <-
   sum(data_HIGH$high.rando == 2 &
         data_HIGH$high.AGE > 68 & data_HIGH$high.dc28 == 1,
       na.rm = T)
 
 
 pi11_HIGH_age_3subset <- r11_age / n11_age
 pi01_HIGH_age_3subset <- r10_age / n10_age
 pi12_HIGH_age_3subset <- r12_age / n12_age
 pi02_HIGH_age_3subset <- r02_age / n02_age
 pi13_HIGH_age_3subset <- r13_age / n13_age
 pi03_HIGH_age_3subset <- r03_age / n03_age
 pa_HIGH_age_sub1 <- mean(data_HIGH$high.AGE <= 58)
 pa_HIGH_age_sub2 <- mean(data_HIGH$high.AGE > 58 & data_HIGH$high.AGE<=68)
 pa_HIGH_age_sub3 <- mean(data_HIGH$high.AGE > 68)
 
 ### functions subset 3 age #####
 simulate_binary_bayesian_trial_3_subsets_HIGH_AGE <-
   function(pi11 ,
            pi01 ,
            pi12 ,
            pi02 ,
            pi13 ,
            pi03 ,
            pa1 ,
            pa2 ,
            pa3 ,
            bayesian.model ,
            stopping_rules.IA ,
            stopping_rules.TA ,
            lambda ,
            gamma ,
            epsilon,
            C1,
            C2) {
     mat.sub1 <- create_matrix_3_subset(Nsimu)
     #for (i in 1:Nsimu) {
     First_result <- NULL
     Second_result <- NULL
     Third_result <- NULL
     Final_result <- NULL
     
     Data_HIGH_IA1 <- filter(data_HIGH, high.Datrand <= IA1)
     subset <- NULL
     subset <- ifelse(Data_HIGH_IA1$high.AGE <= 58, 1, 
                      ifelse(Data_HIGH_IA1$high.AGE>58 &Data_HIGH_IA1$high.AGE<=68,2,
                             ifelse(Data_HIGH_IA1$high.AGE>68,3,0)))
     group <- NULL
     group <- ifelse(Data_HIGH_IA1$high.rando == 1, 1, 0)
     age_11_IA1<- filter(Data_HIGH_IA1, high.rando == 1 &
                            high.AGE <= 58)
     age_01_IA1 <-
       filter(Data_HIGH_IA1, high.rando == 2 & high.AGE <= 58)
     age_12_IA1 <-
       filter(Data_HIGH_IA1, high.rando == 1 & Data_HIGH_IA1$high.AGE > 58 & Data_HIGH_IA1$high.AGE<=68)
     age_02_IA1 <- filter(Data_HIGH_IA1, high.rando == 2 &
                            Data_HIGH_IA1$high.AGE > 58 & Data_HIGH_IA1$high.AGE<=68)
     age_13_IA1 <-
       filter(Data_HIGH_IA1, high.rando == 1 & high.AGE > 68)
     age_03_IA1 <- filter(Data_HIGH_IA1, high.rando == 2 &
                            high.AGE > 68)
     r11_IA1 <- age_11_IA1$high.dc28
     r01_IA1 <- age_01_IA1$high.dc28
     r12_IA1 <- age_12_IA1$high.dc28
     r02_IA1 <- age_02_IA1$high.dc28
     r13_IA1 <- age_13_IA1$high.dc28
     r03_IA1 <- age_03_IA1$high.dc28
     ##### Subset biomarkers
     y <- NULL                #### Response
     y[group == 1 & subset == 1] <- r11_IA1
     y[group == 0 & subset == 1] <- r01_IA1
     y[group == 1 & subset == 2] <- r12_IA1
     y[group == 0 & subset == 2] <- r02_IA1
     y[group == 1 & subset == 3] <- r13_IA1
     y[group == 0 & subset == 3] <- r03_IA1
     Data_IA1 <- data.frame(group, subset, y)
     
     inits_IA1 <-
       create_list_for_inits_values_for_3_subsets(Data_IA1,
                                                  pi11,
                                                  pi01,
                                                  pi12,
                                                  pi02,
                                                  pi13,
                                                  pi03,
                                                  pa1,
                                                  pa2,
                                                  pa3,
                                                  lambda,
                                                  C1,
                                                  C2)
     
     bayesian_jags_fonction_IA1 <-
       jags(
         data = inits_IA1,
         model = bayesian.model,
         n.iter = 30000,
         n.burnin = 20000,
         n.thin = 5,
         n.chains = 3,
         parameters.to.save = c(
           "p1.estim",
           "p0.estim",
           "p11",
           "p01",
           "p12",
           "p02",
           "p13",
           "p03",
           "RR.subset1",
           "RR.subset2",
           "RR.subset3",
           "RR",
           "sd1",
           "sd2",
           "sd3",
           "P1_1",
           "P1_2",
           "P1_3",
           "P5",
           "P6",
           "beta_1",
           "beta_2",
           "beta_3",
           "biais_1",
           "biais_2",
           "biais_3",
           "sum_final"
         ),
         progress.bar = "none",
         jags.seed = "18121995"
       )
     
     
     results_IA1 <-
       results_bayesian_model_3_subsets(bayesian_jags_fonction_IA1, Data_IA1)
     
     First_result <- list(Data_IA1, results_IA1)
     
     Check_interim_decision1 <-
       stopping_rules.IA(
         First_result[[2]][13],
         First_result[[2]][14],
         First_result[[2]][15],
         First_result[[2]][16],
         First_result[[2]][17],
         gamma,
         epsilon
       )
     print(results_IA1)
     print(Check_interim_decision1)
     
     if (Check_interim_decision1[1] == "continue with 3 subsets") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand  <= IA2)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA2$high.AGE <= 58, 1, 
                        ifelse(Data_HIGH_IA2$high.AGE>58 &Data_HIGH_IA2$high.AGE<=68,2,
                               ifelse(Data_HIGH_IA2$high.AGE>68,3,0)))
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_11_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE <= 58)
       age_12_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 1 & Data_HIGH_IA2$high.AGE > 58 & Data_HIGH_IA2$high.AGE<=68)
       age_02_IA2 <- filter(Data_HIGH_IA2, high.rando == 2 &
                              Data_HIGH_IA2$high.AGE > 58 & Data_HIGH_IA2$high.AGE<=68)
       age_13_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 1 & high.AGE > 68)
       age_03_IA2 <- filter(Data_HIGH_IA2, high.rando == 2 &
                              high.AGE > 68)
       r11_IA2 <- age_11_IA2$high.dc28
       r01_IA2 <- age_01_IA2$high.dc28
       r12_IA2 <- age_12_IA2$high.dc28
       r02_IA2 <- age_02_IA2$high.dc28
       r13_IA2 <- age_13_IA2$high.dc28
       r03_IA2 <- age_03_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA2
       y[group == 0 & subset == 1] <- r01_IA2
       y[group == 1 & subset == 2] <- r12_IA2
       y[group == 0 & subset == 2] <- r02_IA2
       y[group == 1 & subset == 3] <- r13_IA2
       y[group == 0 & subset == 3] <- r03_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
       
       
     } else if (Check_interim_decision1[1] == "efficacy in the subset 1 & 2") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand > IA1 &
                                 high.Datrand <= IA2,data_HIGH$high.AGE<=68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA2$high.AGE <= 58, 1, 
                        ifelse(Data_HIGH_IA2$high.AGE>58 &Data_HIGH_IA2$high.AGE<=68,2,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_11_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE <= 58)
       age_12_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 1 & Data_HIGH_IA2$high.AGE > 58 & Data_HIGH_IA2$high.AGE<=68)
       age_02_IA2 <- filter(Data_HIGH_IA2, high.rando == 2 &
                              Data_HIGH_IA2$high.AGE > 58 & Data_HIGH_IA2$high.AGE<=68)
       r11_IA2 <- age_11_IA2$high.dc28
       r01_IA2 <- age_01_IA2$high.dc28
       r12_IA2 <- age_12_IA2$high.dc28
       r02_IA2 <- age_02_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA2
       y[group == 0 & subset == 1] <- r01_IA2
       y[group == 1 & subset == 2] <- r12_IA2
       y[group == 0 & subset == 2] <- r02_IA2
       Data_IA2 <- data.frame(group, subset, y)

       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
       
       
     } else if (Check_interim_decision1[1] == "efficacy in the subset 1 & 3") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand > IA1 &
                                 high.Datrand <= IA2,data_HIGH$high.AGE<=58 &data_HIGH$high.AGE>68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA2$high.AGE <= 58, 1, 
                       ifelse(Data_HIGH_IA2$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_11_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE <= 58)
       age_13_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 1 & Data_HIGH_IA2$high.AGE > 68 )
       age_03_IA2 <- filter(Data_HIGH_IA2, high.rando == 2 &
                              Data_HIGH_IA2$high.AGE > 68 )
       r11_IA2 <- age_11_IA2$high.dc28
       r01_IA2 <- age_01_IA2$high.dc28
       r13_IA2 <- age_13_IA2$high.dc28
       r03_IA2 <- age_03_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA2
       y[group == 0 & subset == 1] <- r01_IA2
       y[group == 1 & subset == 3] <- r13_IA2
       y[group == 0 & subset == 3] <- r03_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
       
       
     } else if (Check_interim_decision1[1] == "efficacy in the subset 2 & 3") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand > IA1 &
                                 high.Datrand <= IA2,data_HIGH$high.AGE>58)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA2$high.AGE> 58& Data_HIGH_IA2$high.AGE<= 68, 2, 
                       ifelse(Data_HIGH_IA2$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_12_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE > 58& high.AGE<=68)
       age_02_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 &  high.AGE > 58 & high.AGE<=68)
       age_13_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 1 & Data_HIGH_IA2$high.AGE > 68 )
       age_03_IA2 <- filter(Data_HIGH_IA2, high.rando == 2 &
                              Data_HIGH_IA2$high.AGE > 68 )
       r12_IA2 <- age_12_IA2$high.dc28
       r02_IA2 <- age_02_IA2$high.dc28
       r13_IA2 <- age_13_IA2$high.dc28
       r03_IA2 <- age_03_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA2
       y[group == 0 & subset == 2] <- r02_IA2
       y[group == 1 & subset == 3] <- r13_IA2
       y[group == 0 & subset == 3] <- r03_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
       
       
     } else if (Check_interim_decision1[1] == "efficacy in the subset 1") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand  > IA1& high.Datrand  <= IA2, high.AGE<=58)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA2$high.AGE <= 58, 1, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_11_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE <= 58)
       r11_IA2 <- age_11_IA2$high.dc28
       r01_IA2 <- age_01_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA2
       y[group == 0 & subset == 1] <- r01_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
       
     } else if (Check_interim_decision1[1] == "efficacy in the subset 2") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand  > IA1& high.Datrand  <= IA2, high.AGE>58& high.AGE<=68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA2$high.AGE > 58 &Data_HIGH_IA2$high.AGE <= 68, 2, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_12_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE > 58 & high.AGE <=68)
       age_02_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE > 58 & high.AGE<=68)
       r12_IA2 <- age_12_IA2$high.dc28
       r02_IA2 <- age_02_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA2
       y[group == 0 & subset == 2] <- r02_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
     } else if (Check_interim_decision1[1] == "efficacy in the subset 3") {
       Data_HIGH_IA2 <- filter(data_HIGH, high.Datrand  > IA1& high.Datrand  <= IA2, high.AGE>68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA2$high.AGE > 68 , 3, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA2$high.rando == 1, 1, 0)
       age_13_IA2<- filter(Data_HIGH_IA2, high.rando == 1 &
                             high.AGE > 68 )
       age_03_IA2 <-
         filter(Data_HIGH_IA2, high.rando == 2 & high.AGE > 68 )
       r13_IA2 <- age_13_IA2$high.dc28
       r03_IA2 <- age_03_IA2$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 3] <- r13_IA2
       y[group == 0 & subset == 3] <- r03_IA2
       Data_IA2 <- data.frame(group, subset, y)
       
       Data_IA2 <- rbind(First_result[[1]], Data_IA2)
       inits_IA2 <-
         create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA2 <-
         jags(
           data = inits_IA2,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA2 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA2, Data_IA2)
       
       Second_result <- list(Data_IA2, results_IA2)
     }
     
     Check_interim_decision2 <-
       stopping_rules.TA(
         Check_interim_decision1[1],
         Second_result[[2]][13],
         Second_result[[2]][14],
         Second_result[[2]][15],
         Second_result[[2]][16],
         Second_result[[2]][17],
         gamma,
         epsilon,
         Check_interim_decision1[14]
       )
     print(results_IA2)
     print(Check_interim_decision2)
     if (Check_interim_decision2[1] == "continue with 3 subsets") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand  <= IA3)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA3$high.AGE <= 58, 1, 
                        ifelse(Data_HIGH_IA3$high.AGE>58 &Data_HIGH_IA3$high.AGE<=68,2,
                               ifelse(Data_HIGH_IA3$high.AGE>68,3,0)))
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_11_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE <= 58)
       age_12_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 1 & Data_HIGH_IA3$high.AGE > 58 & Data_HIGH_IA3$high.AGE<=68)
       age_02_IA3 <- filter(Data_HIGH_IA3, high.rando == 2 &
                              Data_HIGH_IA3$high.AGE > 58 & Data_HIGH_IA3$high.AGE<=68)
       age_13_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 1 & high.AGE > 68)
       age_03_IA3 <- filter(Data_HIGH_IA3, high.rando == 2 &
                              high.AGE > 68)
       r11_IA3 <- age_11_IA3$high.dc28
       r01_IA3 <- age_01_IA3$high.dc28
       r12_IA3 <- age_12_IA3$high.dc28
       r02_IA3 <- age_02_IA3$high.dc28
       r13_IA3 <- age_13_IA3$high.dc28
       r03_IA3 <- age_03_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA3
       y[group == 0 & subset == 1] <- r01_IA3
       y[group == 1 & subset == 2] <- r12_IA3
       y[group == 0 & subset == 2] <- r02_IA3
       y[group == 1 & subset == 3] <- r13_IA3
       y[group == 0 & subset == 3] <- r03_IA3
       
       Data_IA3 <- data.frame(group, subset, y)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
       
       
     } else if (Check_interim_decision2[1] == "efficacy in the subset 1 & 2") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand > IA2 &
                                 high.Datrand <= IA3,data_HIGH$high.AGE<=68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA3$high.AGE <= 58, 1, 
                       ifelse(Data_HIGH_IA3$high.AGE>58 &Data_HIGH_IA3$high.AGE<=68,2,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_11_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE <= 58)
       age_12_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 1 & Data_HIGH_IA3$high.AGE > 58 & Data_HIGH_IA3$high.AGE<=68)
       age_02_IA3 <- filter(Data_HIGH_IA3, high.rando == 2 &
                              Data_HIGH_IA3$high.AGE > 58 & Data_HIGH_IA3$high.AGE<=68)
       r11_IA3 <- age_11_IA3$high.dc28
       r01_IA3 <- age_01_IA3$high.dc28
       r12_IA3 <- age_12_IA3$high.dc28
       r02_IA3 <- age_02_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA3
       y[group == 0 & subset == 1] <- r01_IA3
       y[group == 1 & subset == 2] <- r12_IA3
       y[group == 0 & subset == 2] <- r02_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
       
       
     } else if (Check_interim_decision2[1] == "efficacy in the subset 1 & 3") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand > IA2 &
                                 high.Datrand <= IA3,data_HIGH$high.AGE<=58 &data_HIGH$high.AGE>68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA3$high.AGE <= 58, 1, 
                       ifelse(Data_HIGH_IA3$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_11_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE <= 58)
       age_13_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 1 & Data_HIGH_IA3$high.AGE > 68 )
       age_03_IA3 <- filter(Data_HIGH_IA3, high.rando == 2 &
                              Data_HIGH_IA3$high.AGE > 68 )
       r11_IA3 <- age_11_IA3$high.dc28
       r01_IA3 <- age_01_IA3$high.dc28
       r13_IA3 <- age_13_IA3$high.dc28
       r03_IA3 <- age_03_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA3
       y[group == 0 & subset == 1] <- r01_IA3
       y[group == 1 & subset == 3] <- r13_IA3
       y[group == 0 & subset == 3] <- r03_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
       
       
     } else if (Check_interim_decision2[1] == "efficacy in the subset 2 & 3") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand > IA2 &
                                 high.Datrand <= IA3,data_HIGH$high.AGE>58)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA3$high.AGE> 58& Data_HIGH_IA3$high.AGE<= 68, 2, 
                       ifelse(Data_HIGH_IA3$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_12_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE > 58& high.AGE<=68)
       age_02_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 &  high.AGE > 58 & high.AGE<=68)
       age_13_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 1 & Data_HIGH_IA3$high.AGE > 68 )
       age_03_IA3 <- filter(Data_HIGH_IA3, high.rando == 2 &
                              Data_HIGH_IA3$high.AGE > 68 )
       r12_IA3 <- age_12_IA3$high.dc28
       r02_IA3 <- age_02_IA3$high.dc28
       r13_IA3 <- age_13_IA3$high.dc28
       r03_IA3 <- age_03_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA3
       y[group == 0 & subset == 2] <- r02_IA3
       y[group == 1 & subset == 3] <- r13_IA3
       y[group == 0 & subset == 3] <- r03_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
       
       
     } else if (Check_interim_decision2[1] == "efficacy in the subset 1") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand  > IA2& high.Datrand  <= IA3, high.AGE<=58)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA3$high.AGE <= 58, 1, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_11_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE <= 58)
       r11_IA3 <- age_11_IA3$high.dc28
       r01_IA3 <- age_01_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA3
       y[group == 0 & subset == 1] <- r01_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
       
     } else if (Check_interim_decision2[1] == "efficacy in the subset 2") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand  > IA2& high.Datrand  <= IA3, high.AGE>58& high.AGE<=68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA3$high.AGE > 58 &Data_HIGH_IA3$high.AGE <= 68, 2, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_12_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE > 58 & high.AGE <=68)
       age_02_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE > 58 & high.AGE<=68)
       r12_IA3 <- age_12_IA3$high.dc28
       r02_IA3 <- age_02_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA3
       y[group == 0 & subset == 2] <- r02_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
     } else if (Check_interim_decision2[1] == "efficacy in the subset 3") {
       Data_HIGH_IA3 <- filter(data_HIGH, high.Datrand  > IA2& high.Datrand  <= IA3, high.AGE>68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA3$high.AGE > 68 , 3, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA3$high.rando == 1, 1, 0)
       age_13_IA3<- filter(Data_HIGH_IA3, high.rando == 1 &
                             high.AGE > 68 )
       age_03_IA3 <-
         filter(Data_HIGH_IA3, high.rando == 2 & high.AGE > 68 )
       r13_IA3 <- age_13_IA3$high.dc28
       r03_IA3 <- age_03_IA3$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 3] <- r13_IA3
       y[group == 0 & subset == 3] <- r03_IA3
       Data_IA3 <- data.frame(group, subset, y)
       
       Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
       inits_IA3 <-
         create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA3 <-
         jags(
           data = inits_IA3,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA3 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA3, Data_IA3)
       
       Third_result <- list(Data_IA3, results_IA3)
     }
     
     Check_interim_decision3 <-
       stopping_rules.TA(
         Check_interim_decision2[1],
         Third_result[[2]][13],
         Third_result[[2]][14],
         Third_result[[2]][15],
         Third_result[[2]][16],
         Third_result[[2]][17],
         gamma,
         epsilon,
         Check_interim_decision2[14]
       )
     print(results_IA3)
     print(Check_interim_decision3)
     
     if (Check_interim_decision3[1] == "continue with 3 subsets") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand  <= TA)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA4$high.AGE <= 58, 1, 
                        ifelse(Data_HIGH_IA4$high.AGE>58 &Data_HIGH_IA4$high.AGE<=68,2,
                               ifelse(Data_HIGH_IA4$high.AGE>68,3,0)))
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_11_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE <= 58)
       age_12_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 1 & Data_HIGH_IA4$high.AGE > 58 & Data_HIGH_IA4$high.AGE<=68)
       age_02_IA4 <- filter(Data_HIGH_IA4, high.rando == 2 &
                              Data_HIGH_IA4$high.AGE > 58 & Data_HIGH_IA4$high.AGE<=68)
       age_13_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 1 & high.AGE > 68)
       age_03_IA4 <- filter(Data_HIGH_IA4, high.rando == 2 &
                              high.AGE > 68)
       r11_IA4 <- age_11_IA4$high.dc28
       r01_IA4 <- age_01_IA4$high.dc28
       r12_IA4 <- age_12_IA4$high.dc28
       r02_IA4 <- age_02_IA4$high.dc28
       r13_IA4 <- age_13_IA4$high.dc28
       r03_IA4 <- age_03_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA4
       y[group == 0 & subset == 1] <- r01_IA4
       y[group == 1 & subset == 2] <- r12_IA4
       y[group == 0 & subset == 2] <- r02_IA4
       y[group == 1 & subset == 3] <- r13_IA4
       y[group == 0 & subset == 3] <- r03_IA4
       
       Data_IA4 <- data.frame(group, subset, y)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       
       
     } else if (Check_interim_decision3[1] == "efficacy in the subset 1 & 2") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand > IA3 &
                                 high.Datrand <= TA,data_HIGH$high.AGE<=68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA4$high.AGE <= 58, 1, 
                       ifelse(Data_HIGH_IA4$high.AGE>58 &Data_HIGH_IA4$high.AGE<=68,2,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_11_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE <= 58)
       age_12_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 1 & Data_HIGH_IA4$high.AGE > 58 & Data_HIGH_IA4$high.AGE<=68)
       age_02_IA4 <- filter(Data_HIGH_IA4, high.rando == 2 &
                              Data_HIGH_IA4$high.AGE > 58 & Data_HIGH_IA4$high.AGE<=68)
       r11_IA4 <- age_11_IA4$high.dc28
       r01_IA4 <- age_01_IA4$high.dc28
       r12_IA4 <- age_12_IA4$high.dc28
       r02_IA4 <- age_02_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA4
       y[group == 0 & subset == 1] <- r01_IA4
       y[group == 1 & subset == 2] <- r12_IA4
       y[group == 0 & subset == 2] <- r02_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       
       
     } else if (Check_interim_decision3[1] == "efficacy in the subset 1 & 3") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand > IA3 &
                                 high.Datrand <= TA,data_HIGH$high.AGE<=58 &data_HIGH$high.AGE>68)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA4$high.AGE <= 58, 1, 
                       ifelse(Data_HIGH_IA4$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_11_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE <= 58)
       age_13_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 1 & Data_HIGH_IA4$high.AGE > 68 )
       age_03_IA4 <- filter(Data_HIGH_IA4, high.rando == 2 &
                              Data_HIGH_IA4$high.AGE > 68 )
       r11_IA4 <- age_11_IA4$high.dc28
       r01_IA4 <- age_01_IA4$high.dc28
       r13_IA4 <- age_13_IA4$high.dc28
       r03_IA4 <- age_03_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA4
       y[group == 0 & subset == 1] <- r01_IA4
       y[group == 1 & subset == 3] <- r13_IA4
       y[group == 0 & subset == 3] <- r03_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       
       
     } else if (Check_interim_decision3[1] == "efficacy in the subset 2 & 3") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand > IA3 &
                                 high.Datrand <= TA,data_HIGH$high.AGE>58)
       subset <- NULL
       subset <-ifelse(Data_HIGH_IA4$high.AGE> 58& Data_HIGH_IA4$high.AGE<= 68, 2, 
                       ifelse(Data_HIGH_IA4$high.AGE>68 ,3,0))
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_12_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE > 58& high.AGE<=68)
       age_02_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 &  high.AGE > 58 & high.AGE<=68)
       age_13_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 1 & Data_HIGH_IA4$high.AGE > 68 )
       age_03_IA4 <- filter(Data_HIGH_IA4, high.rando == 2 &
                              Data_HIGH_IA4$high.AGE > 68 )
       r12_IA4 <- age_12_IA4$high.dc28
       r02_IA4 <- age_02_IA4$high.dc28
       r13_IA4 <- age_13_IA4$high.dc28
       r03_IA4 <- age_03_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA4
       y[group == 0 & subset == 2] <- r02_IA4
       y[group == 1 & subset == 3] <- r13_IA4
       y[group == 0 & subset == 3] <- r03_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       
       
     } else if (Check_interim_decision3[1] == "efficacy in the subset 1") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand  > IA3& high.Datrand  <= TA, high.AGE<=58)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA4$high.AGE <= 58, 1, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_11_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE <= 58)
       age_01_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE <= 58)
       r11_IA4 <- age_11_IA4$high.dc28
       r01_IA4 <- age_01_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 1] <- r11_IA4
       y[group == 0 & subset == 1] <- r01_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       
     } else if (Check_interim_decision3[1] == "efficacy in the subset 2") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand  > IA3& high.Datrand  <= TA, high.AGE>58& high.AGE<=68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA4$high.AGE > 58 &Data_HIGH_IA4$high.AGE <= 68, 2, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_12_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE > 58 & high.AGE <=68)
       age_02_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE > 58 & high.AGE<=68)
       r12_IA4 <- age_12_IA4$high.dc28
       r02_IA4 <- age_02_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 2] <- r12_IA4
       y[group == 0 & subset == 2] <- r02_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
     } else if (Check_interim_decision3[1] == "efficacy in the subset 3") {
       Data_HIGH_IA4 <- filter(data_HIGH, high.Datrand  > IA3& high.Datrand  <= TA, high.AGE>68)
       subset <- NULL
       subset <- ifelse(Data_HIGH_IA4$high.AGE > 68 , 3, 0)
       group <- NULL
       group <- ifelse(Data_HIGH_IA4$high.rando == 1, 1, 0)
       age_13_IA4<- filter(Data_HIGH_IA4, high.rando == 1 &
                             high.AGE > 68 )
       age_03_IA4 <-
         filter(Data_HIGH_IA4, high.rando == 2 & high.AGE > 68 )
       r13_IA4 <- age_13_IA4$high.dc28
       r03_IA4 <- age_03_IA4$high.dc28
       ##### Subset biomarkers
       y <- NULL                #### Response
       y[group == 1 & subset == 3] <- r13_IA4
       y[group == 0 & subset == 3] <- r03_IA4
       Data_IA4 <- data.frame(group, subset, y)
       
       Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
       inits_IA4 <-
         create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                    pi11,
                                                    pi01,
                                                    pi12,
                                                    pi02,
                                                    pi13,
                                                    pi03,
                                                    pa1,
                                                    pa2,
                                                    pa3,
                                                    lambda,
                                                    C1,
                                                    C2)
       
       bayesian_jags_fonction_IA4 <-
         jags(
           data = inits_IA4,
           model = bayesian.model,
           n.iter = 30000,
           n.burnin = 20000,
           n.thin = 5,
           n.chains = 3,
           parameters.to.save = c(
             "p1.estim",
             "p0.estim",
             "p11",
             "p01",
             "p12",
             "p02",
             "p13",
             "p03",
             "RR.subset1",
             "RR.subset2",
             "RR.subset3",
             "RR",
             "sd1",
             "sd2",
             "sd3",
             "P1_1",
             "P1_2",
             "P1_3",
             "P5",
             "P6",
             "beta_1",
             "beta_2",
             "beta_3",
             "biais_1",
             "biais_2",
             "biais_3",
             "sum_final"
           ),
           progress.bar = "none",
           jags.seed = "18121995"
         )
       
       
       results_IA4 <-
         results_bayesian_model_3_subsets(bayesian_jags_fonction_IA4, Data_IA4)
       
       Final_result <- list(Data_IA4, results_IA4)
       print(result_IA4)
     }
     
     Check_interim_decision4 <-
       stopping_rules.TA(
         Check_interim_decision3[1],
         Final_result[[2]][13],
         Final_result[[2]][14],
         Final_result[[2]][15],
         Final_result[[2]][16],
         Final_result[[2]][17],
         gamma,
         epsilon,
         Check_interim_decision3[14]
       )
     print(results_IA4)
     print(Check_interim_decision4)
     
     mat.sub1 <-
       cbind(
         gamma,
         C1,
         C2,
         q1,
         q2,
         q3,
         pa1,
         pa2,
         pa3,
         as.numeric(Final_result[[2]][1]),
         as.numeric(Final_result[[2]][2]),
         as.numeric(Final_result[[2]][3]),
         as.numeric(Final_result[[2]][4]),
         as.numeric(Final_result[[2]][5]),
         as.numeric(Final_result[[2]][6]),
         as.numeric(Final_result[[2]][7]),
         as.numeric(Final_result[[2]][8]),
         as.numeric(Final_result[[2]][9]),
         as.numeric(Final_result[[2]][10]),
         as.numeric(Final_result[[2]][11]),
         as.numeric(Final_result[[2]][12]),
         as.numeric(Check_interim_decision1[2]),
         as.numeric(Check_interim_decision1[3]),
         as.numeric(Check_interim_decision1[4]),
         as.numeric(Check_interim_decision1[5]),
         as.numeric(Check_interim_decision1[6]),
         as.numeric(Check_interim_decision1[7]),
         as.numeric(Check_interim_decision1[8]),
         as.numeric(Check_interim_decision1[9]),
         as.numeric(Check_interim_decision1[10]),
         as.numeric(Check_interim_decision1[11]),
         as.numeric(Check_interim_decision1[12]),
         as.numeric(Check_interim_decision1[13]),
         as.numeric(Check_interim_decision2[2]),
         as.numeric(Check_interim_decision2[3]),
         as.numeric(Check_interim_decision2[4]),
         as.numeric(Check_interim_decision2[5]),
         as.numeric(Check_interim_decision2[6]),
         as.numeric(Check_interim_decision2[7]),
         as.numeric(Check_interim_decision2[8]),
         as.numeric(Check_interim_decision2[9]),
         as.numeric(Check_interim_decision2[10]),
         as.numeric(Check_interim_decision2[11]),
         as.numeric(Check_interim_decision2[12]),
         as.numeric(Check_interim_decision2[13]),
         as.numeric(Check_interim_decision3[2]),
         as.numeric(Check_interim_decision3[3]),
         as.numeric(Check_interim_decision3[4]),
         as.numeric(Check_interim_decision3[5]),
         as.numeric(Check_interim_decision3[6]),
         as.numeric(Check_interim_decision3[7]),
         as.numeric(Check_interim_decision3[8]),
         as.numeric(Check_interim_decision3[9]),
         as.numeric(Check_interim_decision3[10]),
         as.numeric(Check_interim_decision3[11]),
         as.numeric(Check_interim_decision3[12]),
         as.numeric(Check_interim_decision3[13]),
         as.numeric(Check_interim_decision4[2]),
         as.numeric(Check_interim_decision4[3]),
         as.numeric(Check_interim_decision4[4]),
         as.numeric(Check_interim_decision4[5]),
         as.numeric(Check_interim_decision4[6]),
         as.numeric(Check_interim_decision4[7]),
         as.numeric(Check_interim_decision4[8]),
         as.numeric(Check_interim_decision4[9]),
         as.numeric(Check_interim_decision4[10]),
         as.numeric(Check_interim_decision4[11]),
         as.numeric(Check_interim_decision4[12]),
         as.numeric(Check_interim_decision4[13]),
         as.numeric(First_result[[2]][30]),
         as.numeric(First_result[[2]][31]),
         as.numeric(First_result[[2]][32]),
         as.numeric(Second_result[[2]][30]),
         as.numeric(Second_result[[2]][31]),
         as.numeric(Second_result[[2]][32]),
         as.numeric(Third_result[[2]][30]),
         as.numeric(Third_result[[2]][31]),
         as.numeric(Third_result[[2]][32]),
         as.numeric(Final_result[[2]][30]),
         as.numeric(Final_result[[2]][31]),
         as.numeric(Final_result[[2]][32]),
         as.numeric(Final_result[[2]][24]),
         as.numeric(Final_result[[2]][25]),
         as.numeric(Final_result[[2]][26]),
         as.numeric(Final_result[[2]][27]),
         as.numeric(Final_result[[2]][28]),
         as.numeric(Final_result[[2]][29])
       )
     #}
     
     return(apply(mat.sub1, 2, function(x) {
       mean(x, na.rm = T)
     }))
     
   }
 
 Nsimu=1
Age_3_subset<-simulate_binary_bayesian_trial_3_subsets_HIGH_AGE(pi11_HIGH_age_3subset,pi01_HIGH_age_3subset,
                                                  pi12_HIGH_age_3subset,pi02_HIGH_age_3subset,
                                                  pi13_HIGH_age_3subset,pi03_HIGH_age_3subset,
                                                  pa_HIGH_age_sub1,pa_HIGH_age_sub2,pa_HIGH_age_sub3,
                                                  binary_bayesian_model_3_subset,
                                                  stopping_rules_first_analyse_variation_3_subset,
                                                  stopping_rules_other_analyse_variation_3_subset,
                                                  0.9,0.65,0.10,0.460,0.215)

names(Age_3_subset) <- c(
     "Gamma",
     "C1",
     "C2",
     "q1",
     "q2",
     "q3",
     "pa1",
     "pa2",
     "pa3",
     "RR.global",
     "2.5%",
     "97.5%",
     "RR.subset1",
     "2.5%",
     "97.5%",
     "RR.subset2",
     "2.5%",
     "97.5%",
     "RR.subset3",
     "2.5%",
     "97.5%",
     "continue_1 with subset ",
     "efficacy1_1",
     "efficacy2_1",
     "efficacy3_1",
     "interaction_quali_1",
     "interaction_quanti_1",
     "enrichment1_1",
     "enrichment2_1",
     "enrichment3_1",
     "enrichment_1.2_1",
     "enrichment_1.3_1",
     "enrichment_2.3_1",
     "continue_2 with subset ",
     "efficacy1_2",
     "efficacy2_2",
     "efficacy3_2",
     "interaction_quali_2",
     "interaction_quanti_2",
     "enrichment1_2",
     "enrichment2_2",
     "enrichment3_2",
     "enrichment_1.2_2",
     "enrichment_1.3_2",
     "enrichment_2.3_2",
     "continue_3 with subset ",
     "efficacy1_3",
     "efficacy2_3",
     "efficacy3_3",
     "interaction_quali_3",
     "interaction_quanti_3",
     "enrichment1_3",
     "enrichment2_3",
   "enrichment3_3",
     "enrichment_1.2_3",
     "enrichment_1.3_3",
     "enrichment_2.3_3",
     "continue_4 with subset ",
     "efficacy1_4",
     "efficacy2_4",
     "efficacy3_4",
     "interaction_quali_4",
     "interaction_quanti_4",
     "enrichment1_4",
     "enrichment2_4",
     "enrichment3_4",
     "enrichment_1.2_4",
     "enrichment_1.3_4",
     "enrichment_2.3_4",
     "N_subset1_stage1",
     "N_subset2_stage1",
     "N_subset3_stage1",
     "N_subset1_stage2",
     "N_subset2_stage2",
     "N_subset3_stage2",
     "N_subset1_stage3",
     "N_subset2_stage3",
   "N_subset3_stage3",
     "N_subset1_stage4",
     "N_subset2_stage4",
     "N_subset3_stage4",
     "biais_1",
     "biais_1_sd",
     "biais_2",
     "biais_2_sd",
     "biais_3",
     "biais_3_sd"
  
   )
Age_3subset


