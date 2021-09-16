install.packages("R2jags")
install.packages("dplyr")
install.packages("doSNOW")
library(R2jags)
library(dplyr)
library(parallel)
library(foreach)
library(doSNOW)

### values explication 
#pi1A= death proportion in subset A and group 1
#pi1B= death proportion in subset B and group 1
#pi0A= death proportion in subset A and group 0
#pi0B= death proportion in subset B and group 0
#lambda_1 = clinical relevant threshold for efficacy measure
#eta_1 = clinical relevant threshold for interaction measure 
#gamma_1 = decision threshold for efficacy measure
#espilon_1 = decision threshold fo interaction measure 
#Q1 = balance of randomization
#pa = prevalence of subset A
#RR.subset 1 = relative risk in subset A
#RR.subset 0 = relative risk in subset B
#N_1,N_2,N_3,N_4 = sample size recruitment at each interim and terminal analysis



####functions#####
# Stopping rules function for first analysis with 2 subset and Millen method
stopping_rules_first_analyse_variation_2_subset_Millen <-
  function(P1_0, P1_1, P5_0, P5_1, gamma_1,epsilon_1) {
    nDecision <- "continue with 2 subsets"
    efficacy_1 <- ifelse(P1_1 > gamma_1, 1, 0)
    efficacy_0 <- ifelse(P1_0 > gamma_1, 1, 0)
    interaction_1 <- ifelse(P5_1 > epsilon_1, 1, 0)
    interaction_0 <- ifelse(P5_0 > epsilon_1, 1, 0)
    enrichment_1 <- ifelse(efficacy_1 == 1 & interaction_0 == 1, 1, 0)
    enrichment_0 <- ifelse(efficacy_0 == 1 & interaction_1 == 1, 1, 0)
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

# Stopping rules function for 2,3 and terminal analysis with 2 subset and Millen method
stopping_rules_other_analyse_variation_2_subset_Millen <-
  function(nDecision,
           P1_0,
           P1_1,
           P5_0,
           P5_1,
           gamma_1,epsilon_1) {
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
      enrichment_1 <- ifelse(efficacy_1 == 1 & interaction_0 == 1, 1, 0)
      enrichment_0 <- ifelse(efficacy_0 == 1 & interaction_1 == 1, 1, 0)
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
#function  for create matrix for data simulation
create_matrix_2.0 <- function(nrow) {
  mat.2 <-
    matrix(nrow = nrow,
           ncol = 55,
           dimnames = list(
             NULL,
             c(
               "q1",
               "pa",
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
               "N_subset0_stage4",
               "biais_0",
               "biais_0_sd",
               "biais_1",
               "biais_1_sd",
               "biais_int_0",
               "biais_int_0_sd",
               "biais_int_1",
               "biais_int_1_sd"
             )
           ))
  
}# create matrix for simulation

## functions for data generation simulation
create_data_frame_for_continue_with_2_subsets <- function(q1,
                                                          N_patients,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa) {
  subset <- rbinom(N_patients, 1, pa) ### subsets
  q0 <- (0.5 - (q1 * pa)) / (1 - pa) ### randomization
  group <- NULL
  group[subset == 1] <- rbinom(sum(subset == 1), 1, q1)
  group[subset == 0] <-
    rbinom(sum(subset == 0), 1, q0) #### Group traitement vs Placebo
  ##### Subset biomarkers
  y <- NULL                 #### Response
  y[group == 1 &
      subset == 1] <-
    rbinom(sum(group == 1 & subset == 1), 1, pi1A)
  y[group == 1 &
      subset == 0] <-
    rbinom(sum(group == 1 & subset == 0), 1, pi1B)
  y[group == 0 &
      subset == 1] <-
    rbinom(sum(group == 0 & subset == 1), 1, pi0A)
  y[group == 0 &
      subset == 0] <-
    rbinom(sum(group == 0 & subset == 0), 1, pi0B)
  
  return(data.frame(group, subset, y))
} # create data if decision rules aren't met

create_data_frame_if_efficacy_subset0 <- function(q1,
                                                  N_patients,
                                                  pi1A,
                                                  pi1B,
                                                  pi0A,
                                                  pi0B,
                                                  pa) {
  subset <- c(rep(0, N_patients))
  q0 <- (0.5 - (q1 * pa)) / (1 - pa)
  group <- NULL
  group[subset == 0] <-
    rbinom(sum(subset == 0), 1, q0) #### Group traitement vs Placebo
  ##### Subset biomarkers
  y <- NULL                 #### Response
  y[group == 1 &
      subset == 0] <-
    rbinom(sum(group == 1 & subset == 0), 1, pi1B)
  y[group == 0 &
      subset == 0] <-
    rbinom(sum(group == 0 & subset == 0), 1, pi0B)
  return(data.frame(group, subset, y))
} # create data if subset 0 is effective

create_data_frame_if_efficacy_subset1 <- function(q1,
                                                  N_patients,
                                                  pi1A,
                                                  pi1B,
                                                  pi0A,
                                                  pi0B,
                                                  pa) {
  subset <- c(rep(1, N_patients))
  q0 <- (0.5 - (q1 * pa)) / (1 - pa)
  group <- NULL
  group[subset == 1] <-
    rbinom(sum(subset == 1), 1, q1) #### Group traitement vs Placebo
  ##### Subset biomarkers
  y <- NULL                 #### Response
  y[group == 1 &
      subset == 1] <-
    rbinom(sum(group == 1 & subset == 1), 1, pi1A)
  y[group == 0 &
      subset == 1] <-
    rbinom(sum(group == 0 & subset == 1), 1, pi0A)
  return(data.frame(group, subset, y))
} # create data if subset 1 is effective


#bayesian model with 2 subset and Millen method for JAGS 
binary_bayesian_model_2_subset_Millen <- function() {
  nb11 ~ dbin(p11, n11)
  nb10 ~ dbin(p10, n10)
  nb01 ~ dbin(p01, n01)
  nb00 ~ dbin(p00, n00)
  
  p10 ~ dbeta(1, 1)
  p00 ~ dbeta(1, 1)
  p01 ~ dbeta(1, 1)
  p11 ~ dbeta(1, 1)
  
  
  p1.estim<-(p11*pa)+(p10*(1-pa)) ### probability in treatment group
  p0.estim<-(p01*pa)+(p00*(1-pa))  ### probability in control group
  RR<-p1.estim/p0.estim 
  RR.subset1 <- p11 / p01 # effect of treatment in subset 1
  RR.subset0 <- p10 / p00 # effect of treatment in subset 0
  interaction_0 <- RR.subset0 / RR.subset1
  interaction_1 <- RR.subset1 / RR.subset0
  P1_1 <- (RR.subset1 < lambda_1)
  P1_0 <- (RR.subset0 < lambda_1)
  P5_0 <- (interaction_0 > eta_1)
  P5_1 <- (interaction_1 > eta_1)
  biais_0 <- RR.subset0 - (pi1B / pi0B)
  biais_1 <- RR.subset1 - (pi1A / pi0A)
  biais_int_0 <-
    interaction_0 - ((pi1B / pi0B) / (pi1A / pi0A))
  biais_int_1 <-
    interaction_1 - ((pi1A / pi0A) / (pi1B / pi0B))
  
} # bayesian model set 1 for jags function

### function for initialization values for jags
create_list_for_inits_values <-
  function(Data, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1) {
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
      pi1A = pi1A,
      pi1B = pi1B,
      pi0A = pi0A,
      pi0B = pi0B,pa=pa,lambda_1=lambda_1,eta_1=eta_1
    )
  }

# function for JAGS results
results_bayesian_model_2_subset_Millen <-
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
    biais_0 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_0
    biais_0_sd <- bayesian_jags_fonction_IA1$BUGSoutput$summary[8, 2]
    biais_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_1
    biais_1_sd <- bayesian_jags_fonction_IA1$BUGSoutput$summary[9, 2]
    biais_int_0 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_int_0
    biais_int_0_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[10, 2]
    biais_int_1 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_int_1
    biais_int_1_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[11, 2]
    N_subset1_stage1 <- sum(Data$subset == 1)
    N_subset0_stage1 <- sum(Data$subset == 0)
    
    results <-
      cbind(
        RR,RR_2.5,RR_97.5,
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
        biais_0,
        biais_0_sd,
        biais_1,
        biais_1_sd,
        biais_int_0,
        biais_int_0_sd,
        biais_int_1,
        biais_int_1_sd,
        N_subset1_stage1,
        N_subset0_stage1
      )
    results
  }# results from bayesian model


# Simulation
#Simuation with 2 subsets and Millen method
simulate_binary_bayesian_trial_2_subset_Millen <-
  function(q1=0.5,
           N_1=200,
           N_2=200,
           N_3=200,
           N_4=200,
           pi1A=0.3,
           pi1B=0.4,
           pi0A=0.3,
           pi0B=0.4,
           pa=0.5,
           bayesian.model=binary_bayesian_model_2_subset_Millen,
           stopping_rules.IA=stopping_rules_first_analyse_variation_2_subset_Millen,
           stopping_rules.TA=stopping_rules_other_analyse_variation_2_subset_Millen,lambda_1=0.9,eta_1=1.25,gamma_1=0.9,epsilon_1=0.90) {
    mat.sub1 <- create_matrix_2.0(Nsimu)
    #for (i in 1:Nsimu) {
      First_result <- NULL
      Second_result <- NULL
      Third_result <- NULL
      Final_result <- NULL
      
      Data_IA1 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                N_1,
                                                                pi1A,
                                                                pi1B,
                                                                pi0A,
                                                                pi0B,
                                                                pa)
      
      inits_IA1 <-
        create_list_for_inits_values(Data_IA1, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
      
      bayesian_jags_fonction_IA1 <-
        jags(
          data = inits_IA1,
          model = bayesian.model,
          n.iter = 30000,
          n.burnin = 20000,
          n.thin = 5,
          n.chains = 3,
          parameters.to.save = c(
            "p1.estim","p0.estim",
            "p11",
            "p10",
            "p01",
            "p00",
            "RR.subset1",
            "RR.subset0","RR",
            "P1_0",
            "P1_1",
            "P5_0",
            "P5_1",
            "interaction_0",
            "interaction_1",
            "biais_0",
            "biais_1",
            "biais_int_0",
            "biais_int_1"
          ),
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA1 <-
        results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA1, Data_IA1)
      
      First_result <- list(Data_IA1, results_IA1)
      
      Check_interim_decision1 <-
        stopping_rules.IA(
          First_result[[2]][10],
          First_result[[2]][11],
          First_result[[2]][12],
          First_result[[2]][13],gamma_1,epsilon_1
        )
      print(Check_interim_decision1)
      
      if (Check_interim_decision1[1] == "continue with 2 subsets" ){
        Data_IA2 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                  N_2,
                                                                  pi1A,
                                                                  pi1B,
                                                                  pi0A,
                                                                  pi0B,
                                                                  pa)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_iA2 <-
          create_list_for_inits_values(Data_IA2, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA2  <-
          jags(
            data = inits_iA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA2 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA2, Data_IA2)
        
        Second_result <- list(Data_IA2, results_IA2)
        
        
      } else if (Check_interim_decision1[1] == "efficacy in the subset 0") {
        Data_IA2 <- create_data_frame_if_efficacy_subset0(q1,
                                                          N_2,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA2 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA2, Data_IA2)
        
        Second_result <- list(Data_IA2, results_IA2)
        
      } else if (Check_interim_decision1[1] == "efficacy in the subset 1") {
        Data_IA2 <- create_data_frame_if_efficacy_subset1(q1,
                                                          N_2,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA2 <- rbind(First_result[[1]], Data_IA2)
        
        inits_IA2 <-
          create_list_for_inits_values(Data_IA2, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA2 <-
          jags(
            data = inits_IA2,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA2 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA2, Data_IA2)
        
        Second_result <- list(Data_IA2, results_IA2)
      }
      
      Check_interim_decision2 <-
        stopping_rules.TA(
          Check_interim_decision1[1],
          Second_result[[2]][10],
          Second_result[[2]][11],
          Second_result[[2]][12],
          Second_result[[2]][13],gamma_1,epsilon_1
        )
      print(Check_interim_decision2)
      if (
        Check_interim_decision2[1] == "continue with 2 subsets" ) {
        Data_IA3 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                  N_3,
                                                                  pi1A,
                                                                  pi1B,
                                                                  pi0A,
                                                                  pi0B,
                                                                  pa)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA3 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA3, Data_IA3)
        
        Third_result <- list(Data_IA3, results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "efficacy in the subset 0"
      )) {
        Data_IA3 <- create_data_frame_if_efficacy_subset0(q1,
                                                          N_3,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        
        
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA3 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA3, Data_IA3)
        
        Third_result <- list(Data_IA3, results_IA3)
        
      } else if ((
        Check_interim_decision2[1] == "efficacy in the subset 1"
      )) {
        Data_IA3 <- create_data_frame_if_efficacy_subset1(q1,
                                                          N_3,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
        
        inits_IA3 <-
          create_list_for_inits_values(Data_IA3, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA3 <-
          jags(
            data = inits_IA3,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA3 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA3, Data_IA3)
        
        Third_result <- list(Data_IA3, results_IA3)
      }
      
      Check_interim_decision3 <-
        stopping_rules.TA(
          Check_interim_decision2[1],
          Third_result[[2]][10],
          Third_result[[2]][11],
          Third_result[[2]][12],
          Third_result[[2]][13],gamma_1,epsilon_1
        )
      print(Check_interim_decision3)
      if (
        Check_interim_decision3[1] == "continue with 2 subsets" ) {
        Data_IA4 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                  N_4,
                                                                  pi1A,
                                                                  pi1B,
                                                                  pi0A,
                                                                  pi0B,
                                                                  pa)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA4 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA4, Data_IA4)
        
        Final_result <- list(Data_IA4, results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "efficacy in the subset 0"
      )) {
        Data_IA4 <- create_data_frame_if_efficacy_subset0(q1,
                                                          N_4,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA4 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA4, Data_IA4)
        
        Final_result <- list(Data_IA4, results_IA4)
        
      } else if ((
        Check_interim_decision3[1] == "efficacy in the subset 1"
      )) {
        Data_IA4 <- create_data_frame_if_efficacy_subset1(q1,
                                                          N_4,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa)
        Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
        
        inits_IA4 <-
          create_list_for_inits_values(Data_IA4, pi1A, pi1B, pi0A, pi0B,pa,lambda_1,eta_1)
        
        bayesian_jags_fonction_IA4 <-
          jags(
            data = inits_IA4,
            model = bayesian.model,
            n.iter = 30000,
            n.burnin = 20000,
            n.thin = 5,
            n.chains = 3,
            parameters.to.save = c(
              "p1.estim","p0.estim",
              "p11",
              "p10",
              "p01",
              "p00",
              "RR.subset1",
              "RR.subset0","RR",
              "P1_0",
              "P1_1",
              "P5_0",
              "P5_1",
              "interaction_0",
              "interaction_1",
              "biais_0",
              "biais_1",
              "biais_int_0",
              "biais_int_1"
            ),
            progress.bar = "none",
            jags.seed = "18121995"
          )
        
        results_IA4 <-
          results_bayesian_model_2_subset_Millen(bayesian_jags_fonction_IA4, Data_IA4)
        
        Final_result <- list(Data_IA4, results_IA4)
        
        
      }
      
      Check_interim_decision4 <-
        stopping_rules.TA(
          Check_interim_decision3[1],
          Final_result[[2]][10],
          Final_result[[2]][11],
          Final_result[[2]][12],
          Final_result[[2]][13],gamma_1,epsilon_1
        )
      print(Check_interim_decision4)
      
      mat.sub1 <-
        cbind(
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
          as.numeric(First_result[[2]][22]),
          as.numeric(First_result[[2]][23]),
          as.numeric(Second_result[[2]][22]),
          as.numeric(Second_result[[2]][23]),
          as.numeric(Third_result[[2]][22]),
          as.numeric(Third_result[[2]][23]),
          as.numeric(Final_result[[2]][22]),
          as.numeric(Final_result[[2]][23]),
          as.numeric(Final_result[[2]][14]),
          as.numeric(Final_result[[2]][15]),
          as.numeric(Final_result[[2]][16]),
          as.numeric(Final_result[[2]][17]),
          as.numeric(Final_result[[2]][18]),
          as.numeric(Final_result[[2]][19]),
          as.numeric(Final_result[[2]][20]),
          as.numeric(Final_result[[2]][21])
        )
    #}
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }



####simulation with 2 subset and millen interaction method #####
####### values for parallelize the scenarios
Nsimu <- 10000
seeds <- sample(1:200000, Nsimu)
pa_choice <- c(0.2, 0.4, 0.5, 0.6, 0.8)
q1_choice <- c(0.1, 0.3, 0.5, 0.6, 0.9)
cl <- makeCluster(110,outfile="")
registerDoSNOW(cl)


### scenarios with pa( prevalence subset) simulation####
sc1_2subset_variation_pa_Millen <- list()
for (m in seq_along(pa_choice)) {
  pa_choices <- pa_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      0.5,
      200,
      200,
      200,
      200,
      0.4,
      0.4,
      0.4,
      0.4,
      pa_choices,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc1_2subset_variation_pa_Millen[[m]] <- data.frame(res)
}
saveRDS(sc1_2subset_variation_pa_Millen, file = "sc1_2subset_variation_estimate_pa_Millen.rds")
sc1_2subset_variation_pa_Millen <- bind_rows(sc1_2subset_variation_pa_Millen)

sc1_2subset_variation_pa_Millen <-
  sc1_2subset_variation_pa_Millen %>% group_by(pa) %>% summarise_all(mean)
colnames(sc1_2subset_variation_pa_Millen) <- c(
  "pa",
  "q1",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

sc2_2subset_variation_pa_Millen <- list()
for (m in seq_along(pa_choice)) {
  pa_choices <- pa_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      0.5,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa_choices,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc2_2subset_variation_pa_Millen[[m]] <- data.frame(res)
}
saveRDS(sc2_2subset_variation_pa_Millen, file = "sc2_2subset_variation_estimate_pa_Millen.rds")
sc2_2subset_variation_pa_Millen <- bind_rows(sc2_2subset_variation_pa_Millen)

sc2_2subset_variation_pa_Millen <-
  sc2_2subset_variation_pa_Millen %>% group_by(pa) %>% summarise_all(mean)
colnames(sc2_2subset_variation_pa_Millen) <- c(
  "pa",
  "q1",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

sc3_2subset_variation_pa_Millen <- list()
for (m in seq_along(pa_choice)) {
  pa_choices <- pa_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      0.5,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.37,
      0.4,
      pa_choices,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc3_2subset_variation_pa_Millen[[m]] <- data.frame(res)
}
saveRDS(sc3_2subset_variation_pa_Millen, file = "sc3_2subset_variation_estimate_pa_Millen.rds")
sc3_2subset_variation_pa_Millen <- bind_rows(sc3_2subset_variation_pa_Millen)

sc3_2subset_variation_pa_Millen <-
  sc3_2subset_variation_pa_Millen %>% group_by(pa) %>% summarise_all(mean)
colnames(sc3_2subset_variation_pa_Millen) <- c(
  "pa",
  "q1",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)
sc4_2subset_variation_pa_Millen <- list()
for (m in seq_along(pa_choice)) {
  pa_choices <- pa_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      0.5,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa_choices,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc4_2subset_variation_pa_Millen[[m]] <- data.frame(res)
}
saveRDS(sc4_2subset_variation_pa_Millen, file = "sc4_2subset_variation_estimate_pa_Millen.rds")
sc4_2subset_variation_pa_Millen <- bind_rows(sc4_2subset_variation_pa_Millen)

sc4_2subset_variation_pa_Millen <-
  sc4_2subset_variation_pa_Millen %>% group_by(pa) %>% summarise_all(mean)
colnames(sc4_2subset_variation_pa_Millen) <- c(
  "pa",
  "q1",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

#### scenarios with Q1 ( balance of randomization ) simulation ####
sc1_2subset_variation_q1_Millen <- list()
for (m in seq_along(q1_choice)) {
  q1_choices <- q1_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      q1_choices,
      200,
      200,
      200,
      200,
      0.4,
      0.4,
      0.4,
      0.4,
      0.5,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc1_2subset_variation_q1_Millen[[m]] <- data.frame(res)
}
saveRDS(sc1_2subset_variation_q1_Millen, file = "sc1_2subset_variation_estimate_q1_Millen.rds")
sc1_2subset_variation_q1_Millen <- bind_rows(sc1_2subset_variation_q1_Millen)

sc1_2subset_variation_q1_Millen <-
  sc1_2subset_variation_q1_Millen %>% group_by(q1) %>% summarise_all(mean)
colnames(sc1_2subset_variation_q1_Millen) <- c(
  "q1",
  "pa",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

sc2_2subset_variation_q1_Millen <- list()
for (m in seq_along(q1_choice)) {
  q1_choices <- q1_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      q1_choices,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      0.5,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc2_2subset_variation_q1_Millen[[m]] <- data.frame(res)
}
saveRDS(sc2_2subset_variation_q1_Millen, file = "sc2_2subset_variation_estimate_q1_Millen.rds")
sc2_2subset_variation_q1_Millen <- bind_rows(sc2_2subset_variation_q1_Millen)

sc2_2subset_variation_q1_Millen <-
  sc2_2subset_variation_q1_Millen %>% group_by(q1) %>% summarise_all(mean)
colnames(sc2_2subset_variation_q1_Millen) <- c(
  "q1",
  "pa",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

sc3_2subset_variation_q1_Millen <- list()
for (m in seq_along(q1_choice)) {
  q1_choices <- q1_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      q1_choices,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.37,
      0.4,
      0.5,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc3_2subset_variation_q1_Millen[[m]] <- data.frame(res)
}
saveRDS(sc3_2subset_variation_q1_Millen, file = "sc3_2subset_variation_estimate_q1_Millen.rds")
sc3_2subset_variation_q1_Millen <- bind_rows(sc3_2subset_variation_q1_Millen)

sc3_2subset_variation_q1_Millen <-
  sc3_2subset_variation_q1_Millen %>% group_by(q1) %>% summarise_all(mean)
colnames(sc3_2subset_variation_q1_Millen) <- c(
  "q1",
  "pa",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)

sc4_2subset_variation_q1_Millen <- list()
for (m in seq_along(q1_choice)) {
  q1_choices <- q1_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2.0(
      q1_choices,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      0.5,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      1.25,
      0.90,
      0.90
      
      
    )
  }
  sc4_2subset_variation_q1_Millen[[m]] <- data.frame(res)
}
saveRDS(sc4_2subset_variation_q1_Millen, file = "sc4_2subset_variation_estimate_q1_Millen.rds")
sc4_2subset_variation_q1_Millen <- bind_rows(sc4_2subset_variation_q1_Millen)

sc4_2subset_variation_q1_Millen <-
  sc4_2subset_variation_q1_Millen %>% group_by(q1) %>% summarise_all(mean)
colnames(sc4_2subset_variation_q1_Millen) <- c(
  "q1",
  "pa",
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
  "N_subset0_stage4",
  "biais_0",
  "biais_0_sd",
  "biais_1",
  "biais_1_sd",
  "biais_int_0",
  "biais_int_0_sd",
  "biais_int_1",
  "biais_int_1_sd"
)


#save all scenarios for create graphics and tables
saveRDS(sc1_2subset_variation_pa_Millen, file = "sc1_2subset_variation_pa_Millen.rds")
saveRDS(sc2_2subset_variation_pa_Millen, file = "sc2_2subset_variation_pa_Millen.rds")
saveRDS(sc3_2subset_variation_pa_Millen, file = "sc3_2subset_variation_pa_Millen.rds")
saveRDS(sc4_2subset_variation_pa_Millen, file = "sc4_2subset_variation_pa_Millen.rds")

saveRDS(sc1_2subset_variation_q1_Millen, file = "sc1_2subset_variation_q1_Millen.rds")
saveRDS(sc2_2subset_variation_q1_Millen, file = "sc2_2subset_variation_q1_Millen.rds")
saveRDS(sc3_2subset_variation_q1_Millen, file = "sc3_2subset_variation_q1_Millen.rds")
saveRDS(sc4_2subset_variation_q1_Millen, file = "sc4_2subset_variation_q1_Millen.rds")

