install.packages("R2jags")
install.packages("dplyr")
install.packages("doSNOW")
library(R2jags)
library(dplyr)
library(parallel)
library(foreach)
library(doSNOW)

####simulation####
###Load Simulation_method_with_Millen_interaction for all functions

Nsimu=1000
deb <- Sys.time()
cl<-makeCluster(110)
registerDoSNOW(cl)
eta_choices<-seq(1,2,0.05)

sc1_0.90_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.90
    )
  }
  sc1_0.90_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_0.90_epsilon_0.95_gamma,file="sc1_0.90_epsilon_0.95_gamma_estimate.rds")
sc1_0.90_epsilon_0.95_gamma<-bind_rows(sc1_0.90_epsilon_0.95_gamma)

sc1_0.90_epsilon_0.95_gamma<-sc1_0.90_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_0.90_epsilon_0.95_gamma) <- c("eta",
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

sc1_0.90_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.90
    )
  }
  sc1_0.90_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_0.90_epsilon_0.90_gamma,file="sc1_0.90_epsilon_0.90_gamma_estimate.rds")
sc1_0.90_epsilon_0.90_gamma<-bind_rows(sc1_0.90_epsilon_0.90_gamma)

sc1_0.90_epsilon_0.90_gamma<-sc1_0.90_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_0.90_epsilon_0.90_gamma) <- c("eta",
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

sc1_0.90_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.90
    )
  }
  sc1_0.90_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_0.90_epsilon_0.85_gamma,file="sc1_0.90_epsilon_0.85_gamma_estimate.rds")
sc1_0.90_epsilon_0.85_gamma<-bind_rows(sc1_0.90_epsilon_0.85_gamma)

sc1_0.90_epsilon_0.85_gamma<-sc1_0.90_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_0.90_epsilon_0.85_gamma) <- c("eta",
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
sc1_0.90_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.90
    )
  }
  sc1_0.90_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_0.90_epsilon_0.80_gamma,file="sc1_0.90_epsilon_0.80_gamma_estimate.rds")
sc1_0.90_epsilon_0.80_gamma<-bind_rows(sc1_0.90_epsilon_0.80_gamma)

sc1_0.90_epsilon_0.80_gamma<-sc1_0.90_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_0.90_epsilon_0.80_gamma) <- c("eta",
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
sc1_0.90_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.90
    )
  }
  sc1_0.90_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_0.90_epsilon_0.75_gamma,file="sc1_0.90_epsilon_0.75_gamma_estimate.rds")
sc1_0.90_epsilon_0.75_gamma<-bind_rows(sc1_0.90_epsilon_0.75_gamma)

sc1_0.90_epsilon_0.75_gamma<-sc1_0.90_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_0.90_epsilon_0.75_gamma) <- c("eta",
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

sc1_95_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.95
    )
  }
  sc1_95_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_95_epsilon_0.95_gamma,file="sc1_95_epsilon_0.95_gamma_estimate.rds")
sc1_95_epsilon_0.95_gamma<-bind_rows(sc1_95_epsilon_0.95_gamma)

sc1_95_epsilon_0.95_gamma<-sc1_95_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_95_epsilon_0.95_gamma) <- c("eta",
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

sc1_95_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.95
    )
  }
  sc1_95_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_95_epsilon_0.90_gamma,file="sc1_95_epsilon_0.90_gamma_estimate.rds")
sc1_95_epsilon_0.90_gamma<-bind_rows(sc1_95_epsilon_0.90_gamma)

sc1_95_epsilon_0.90_gamma<-sc1_95_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_95_epsilon_0.90_gamma) <- c("eta",
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

sc1_95_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.95
    )
  }
  sc1_95_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_95_epsilon_0.85_gamma,file="sc1_95_epsilon_0.85_gamma_estimate.rds")
sc1_95_epsilon_0.85_gamma<-bind_rows(sc1_95_epsilon_0.85_gamma)

sc1_95_epsilon_0.85_gamma<-sc1_95_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_95_epsilon_0.85_gamma) <- c("eta",
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
sc1_95_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.95
    )
  }
  sc1_95_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_95_epsilon_0.80_gamma,file="sc1_95_epsilon_0.80_gamma_estimate.rds")
sc1_95_epsilon_0.80_gamma<-bind_rows(sc1_95_epsilon_0.80_gamma)

sc1_95_epsilon_0.80_gamma<-sc1_95_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_95_epsilon_0.80_gamma) <- c("eta",
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
sc1_95_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.3,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.95
    )
  }
  sc1_95_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc1_95_epsilon_0.75_gamma,file="sc1_95_epsilon_0.75_gamma_estimate.rds")
sc1_95_epsilon_0.75_gamma<-bind_rows(sc1_95_epsilon_0.75_gamma)

sc1_95_epsilon_0.75_gamma<-sc1_95_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc1_95_epsilon_0.75_gamma) <- c("eta",
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

sc2_0.90_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.90
    )
  }
  sc2_0.90_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_0.90_epsilon_0.95_gamma,file="sc2_0.90_epsilon_0.95_gamma_estimate.rds")
sc2_0.90_epsilon_0.95_gamma<-bind_rows(sc2_0.90_epsilon_0.95_gamma)

sc2_0.90_epsilon_0.95_gamma<-sc2_0.90_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_0.90_epsilon_0.95_gamma) <- c("eta",
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

sc2_0.90_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.90
    )
  }
  sc2_0.90_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_0.90_epsilon_0.90_gamma,file="sc2_0.90_epsilon_0.90_gamma_estimate.rds")
sc2_0.90_epsilon_0.90_gamma<-bind_rows(sc2_0.90_epsilon_0.90_gamma)

sc2_0.90_epsilon_0.90_gamma<-sc2_0.90_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_0.90_epsilon_0.90_gamma) <- c("eta",
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

sc2_0.90_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.90
    )
  }
  sc2_0.90_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_0.90_epsilon_0.85_gamma,file="sc2_0.90_epsilon_0.85_gamma_estimate.rds")
sc2_0.90_epsilon_0.85_gamma<-bind_rows(sc2_0.90_epsilon_0.85_gamma)

sc2_0.90_epsilon_0.85_gamma<-sc2_0.90_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_0.90_epsilon_0.85_gamma) <- c("eta",
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
sc2_0.90_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.90
    )
  }
  sc2_0.90_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_0.90_epsilon_0.80_gamma,file="sc2_0.90_epsilon_0.80_gamma_estimate.rds")
sc2_0.90_epsilon_0.80_gamma<-bind_rows(sc2_0.90_epsilon_0.80_gamma)

sc2_0.90_epsilon_0.80_gamma<-sc2_0.90_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_0.90_epsilon_0.80_gamma) <- c("eta",
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
sc2_0.90_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.90
    )
  }
  sc2_0.90_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_0.90_epsilon_0.75_gamma,file="sc2_0.90_epsilon_0.75_gamma_estimate.rds")
sc2_0.90_epsilon_0.75_gamma<-bind_rows(sc2_0.90_epsilon_0.75_gamma)

sc2_0.90_epsilon_0.75_gamma<-sc2_0.90_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_0.90_epsilon_0.75_gamma) <- c("eta",
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

sc2_95_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.95
    )
  }
  sc2_95_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_95_epsilon_0.95_gamma,file="sc2_95_epsilon_0.95_gamma_estimate.rds")
sc2_95_epsilon_0.95_gamma<-bind_rows(sc2_95_epsilon_0.95_gamma)

sc2_95_epsilon_0.95_gamma<-sc2_95_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_95_epsilon_0.95_gamma) <- c("eta",
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

sc2_95_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.95
    )
  }
  sc2_95_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_95_epsilon_0.90_gamma,file="sc2_95_epsilon_0.90_gamma_estimate.rds")
sc2_95_epsilon_0.90_gamma<-bind_rows(sc2_95_epsilon_0.90_gamma)

sc2_95_epsilon_0.90_gamma<-sc2_95_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_95_epsilon_0.90_gamma) <- c("eta",
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

sc2_95_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.95
    )
  }
  sc2_95_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_95_epsilon_0.85_gamma,file="sc2_95_epsilon_0.85_gamma_estimate.rds")
sc2_95_epsilon_0.85_gamma<-bind_rows(sc2_95_epsilon_0.85_gamma)

sc2_95_epsilon_0.85_gamma<-sc2_95_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_95_epsilon_0.85_gamma) <- c("eta",
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
sc2_95_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.95
    )
  }
  sc2_95_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_95_epsilon_0.80_gamma,file="sc2_95_epsilon_0.80_gamma_estimate.rds")
sc2_95_epsilon_0.80_gamma<-bind_rows(sc2_95_epsilon_0.80_gamma)

sc2_95_epsilon_0.80_gamma<-sc2_95_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_95_epsilon_0.80_gamma) <- c("eta",
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
sc2_95_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.3,
      0.4,
      0.4,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.95
    )
  }
  sc2_95_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc2_95_epsilon_0.75_gamma,file="sc2_95_epsilon_0.75_gamma_estimate.rds")
sc2_95_epsilon_0.75_gamma<-bind_rows(sc2_95_epsilon_0.75_gamma)

sc2_95_epsilon_0.75_gamma<-sc2_95_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc2_95_epsilon_0.75_gamma) <- c("eta",
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


sc4_0.90_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.90
    )
  }
  sc4_0.90_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_0.90_epsilon_0.95_gamma,file="sc4_0.90_epsilon_0.95_gamma_estimate.rds")
sc4_0.90_epsilon_0.95_gamma<-bind_rows(sc4_0.90_epsilon_0.95_gamma)

sc4_0.90_epsilon_0.95_gamma<-sc4_0.90_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_0.90_epsilon_0.95_gamma) <- c("eta",
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

sc4_0.90_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.90
    )
  }
  sc4_0.90_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_0.90_epsilon_0.90_gamma,file="sc4_0.90_epsilon_0.90_gamma_estimate.rds")
sc4_0.90_epsilon_0.90_gamma<-bind_rows(sc4_0.90_epsilon_0.90_gamma)

sc4_0.90_epsilon_0.90_gamma<-sc4_0.90_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_0.90_epsilon_0.90_gamma) <- c("eta",
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

sc4_0.90_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.90
    )
  }
  sc4_0.90_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_0.90_epsilon_0.85_gamma,file="sc4_0.90_epsilon_0.85_gamma_estimate.rds")
sc4_0.90_epsilon_0.85_gamma<-bind_rows(sc4_0.90_epsilon_0.85_gamma)

sc4_0.90_epsilon_0.85_gamma<-sc4_0.90_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_0.90_epsilon_0.85_gamma) <- c("eta",
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
sc4_0.90_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.90
    )
  }
  sc4_0.90_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_0.90_epsilon_0.80_gamma,file="sc4_0.90_epsilon_0.80_gamma_estimate.rds")
sc4_0.90_epsilon_0.80_gamma<-bind_rows(sc4_0.90_epsilon_0.80_gamma)

sc4_0.90_epsilon_0.80_gamma<-sc4_0.90_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_0.90_epsilon_0.80_gamma) <- c("eta",
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
sc4_0.90_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.90
    )
  }
  sc4_0.90_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_0.90_epsilon_0.75_gamma,file="sc4_0.90_epsilon_0.75_gamma_estimate.rds")
sc4_0.90_epsilon_0.75_gamma<-bind_rows(sc4_0.90_epsilon_0.75_gamma)

sc4_0.90_epsilon_0.75_gamma<-sc4_0.90_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_0.90_epsilon_0.75_gamma) <- c("eta",
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

sc4_95_epsilon_0.95_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.90,
      eta_choice,
      0.95,
      0.95
    )
  }
  sc4_95_epsilon_0.95_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_95_epsilon_0.95_gamma,file="sc4_95_epsilon_0.95_gamma_estimate.rds")
sc4_95_epsilon_0.95_gamma<-bind_rows(sc4_95_epsilon_0.95_gamma)

sc4_95_epsilon_0.95_gamma<-sc4_95_epsilon_0.95_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_95_epsilon_0.95_gamma) <- c("eta",
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

sc4_95_epsilon_0.90_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.90,
      0.95
    )
  }
  sc4_95_epsilon_0.90_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_95_epsilon_0.90_gamma,file="sc4_95_epsilon_0.90_gamma_estimate.rds")
sc4_95_epsilon_0.90_gamma<-bind_rows(sc4_95_epsilon_0.90_gamma)

sc4_95_epsilon_0.90_gamma<-sc4_95_epsilon_0.90_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_95_epsilon_0.90_gamma) <- c("eta",
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

sc4_95_epsilon_0.85_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.85,
      0.95
    )
  }
  sc4_95_epsilon_0.85_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_95_epsilon_0.85_gamma,file="sc4_95_epsilon_0.85_gamma_estimate.rds")
sc4_95_epsilon_0.85_gamma<-bind_rows(sc4_95_epsilon_0.85_gamma)

sc4_95_epsilon_0.85_gamma<-sc4_95_epsilon_0.85_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_95_epsilon_0.85_gamma) <- c("eta",
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
sc4_95_epsilon_0.80_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.80,
      0.95
    )
  }
  sc4_95_epsilon_0.80_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_95_epsilon_0.80_gamma,file="sc4_95_epsilon_0.80_gamma_estimate.rds")
sc4_95_epsilon_0.80_gamma<-bind_rows(sc4_95_epsilon_0.80_gamma)

sc4_95_epsilon_0.80_gamma<-sc4_95_epsilon_0.80_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_95_epsilon_0.80_gamma) <- c("eta",
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
sc4_95_epsilon_0.75_gamma<-list()
for(m in seq_along(eta_choices)){
  eta_choice<-eta_choices[m]
  
  res<- foreach(i=1:Nsimu,
                .inorder=F,.combine="rbind",
                .packages=c("dplyr","R2jags")
  ) %dopar% {
    simulate_binary_bayesian_trial_2_subset_Millen(
      q1,
      200,
      200,
      200,
      200,
      0.2,
      0.4,
      0.5,
      0.4,
      pa,
      binary_bayesian_model_2_subset_Millen,
      stopping_rules_first_analyse_variation_2_subset_Millen,
      stopping_rules_other_analyse_variation_2_subset_Millen,
      0.9,
      eta_choice,
      0.75,
      0.95
    )
  }
  sc4_95_epsilon_0.75_gamma[[m]]<- data.frame(res)
}
saveRDS(sc4_95_epsilon_0.75_gamma,file="sc4_95_epsilon_0.75_gamma_estimate.rds")
sc4_95_epsilon_0.75_gamma<-bind_rows(sc4_95_epsilon_0.75_gamma)

sc4_95_epsilon_0.75_gamma<-sc4_95_epsilon_0.75_gamma %>% group_by(eta_1) %>% summarise_all(mean)
colnames(sc4_95_epsilon_0.75_gamma) <- c("eta",
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
saveRDS(sc1_0.90_epsilon_0.95_gamma,file="sc1_0.90_epsilon_0.95_gamma.rds")
saveRDS(sc1_0.90_epsilon_0.90_gamma,file="sc1_0.90_epsilon_0.90_gamma.rds")
saveRDS(sc1_0.90_epsilon_0.85_gamma,file="sc1_0.90_epsilon_0.85_gamma.rds")
saveRDS(sc1_0.90_epsilon_0.80_gamma,file="sc1_0.90_epsilon_0.80_gamma.rds")
saveRDS(sc1_0.90_epsilon_0.75_gamma,file="sc1_0.90_epsilon_0.75_gamma.rds")

saveRDS(sc1_95_epsilon_0.95_gamma,file="sc1_95_epsilon_0.95_gamma.rds")
saveRDS(sc1_95_epsilon_0.90_gamma,file="sc1_95_epsilon_0.90_gamma.rds")
saveRDS(sc1_95_epsilon_0.85_gamma,file="sc1_95_epsilon_0.85_gamma.rds")
saveRDS(sc1_95_epsilon_0.80_gamma,file="sc1_95_epsilon_0.80_gamma.rds")
saveRDS(sc1_95_epsilon_0.75_gamma,file="sc1_95_epsilon_0.75_gamma.rds")

saveRDS(sc2_0.90_epsilon_0.95_gamma,file="sc2_0.90_epsilon_0.95_gamma.rds")
saveRDS(sc2_0.90_epsilon_0.90_gamma,file="sc2_0.90_epsilon_0.90_gamma.rds")
saveRDS(sc2_0.90_epsilon_0.85_gamma,file="sc2_0.90_epsilon_0.85_gamma.rds")
saveRDS(sc2_0.90_epsilon_0.80_gamma,file="sc2_0.90_epsilon_0.80_gamma.rds")
saveRDS(sc2_0.90_epsilon_0.75_gamma,file="sc2_0.90_epsilon_0.75_gamma.rds")

saveRDS(sc2_95_epsilon_0.95_gamma,file="sc2_95_epsilon_0.95_gamma.rds")
saveRDS(sc2_95_epsilon_0.90_gamma,file="sc2_95_epsilon_0.90_gamma.rds")
saveRDS(sc2_95_epsilon_0.85_gamma,file="sc2_95_epsilon_0.85_gamma.rds")
saveRDS(sc2_95_epsilon_0.80_gamma,file="sc2_95_epsilon_0.80_gamma.rds")
saveRDS(sc2_95_epsilon_0.75_gamma,file="sc2_95_epsilon_0.75_gamma.rds")

saveRDS(sc4_0.90_epsilon_0.95_gamma,file="sc4_0.90_epsilon_0.95_gamma.rds")
saveRDS(sc4_0.90_epsilon_0.90_gamma,file="sc4_0.90_epsilon_0.90_gamma.rds")
saveRDS(sc4_0.90_epsilon_0.85_gamma,file="sc4_0.90_epsilon_0.85_gamma.rds")
saveRDS(sc4_0.90_epsilon_0.80_gamma,file="sc4_0.90_epsilon_0.80_gamma.rds")
saveRDS(sc4_0.90_epsilon_0.75_gamma,file="sc4_0.90_epsilon_0.75_gamma.rds")

saveRDS(sc4_95_epsilon_0.95_gamma,file="sc4_95_epsilon_0.95_gamma.rds")
saveRDS(sc4_95_epsilon_0.90_gamma,file="sc4_95_epsilon_0.90_gamma.rds")
saveRDS(sc4_95_epsilon_0.85_gamma,file="sc4_95_epsilon_0.85_gamma.rds")
saveRDS(sc4_95_epsilon_0.80_gamma,file="sc4_95_epsilon_0.80_gamma.rds")
saveRDS(sc4_95_epsilon_0.75_gamma,file="sc4_95_epsilon_0.75_gamma.rds")
fin <- Sys.time() - deb
stopCluster(cl)