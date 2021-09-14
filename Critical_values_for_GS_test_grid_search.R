install.packages("R2jags")
install.packages("dplyr")
install.packages("doSNOW")
library(R2jags)
library(dplyr)
library(parallel)
library(foreach)
library(doSNOW)
library(ggplot2)
library(ggsci)
###load Simulation_method_with_GS_interaction and Simulations_method_with_Millen_interaction
####simulation  for thresholds parameters for gail and Simon  method#####
Nsimu=1000
seeds <- sample(1:200000,Nsimu)
C_var<-seq(0,2,0.005)
cl<-makeCluster(100,outfile="")
registerDoSNOW(cl)

####Scenarios with 3 subset and under the null hypothesis ####
 sc1_3subset_variation_0.95<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]

   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,0.5,0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       1/3,1/3,1/3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,0.9,
       0.05,
       C_vars

     )
   }
   sc1_3subset_variation_0.95[[m]]<- data.frame(res)
 }
 saveRDS(sc1_3subset_variation_0.95,file="sc1_3subset_variation_estimate_0.95.rds")
 sc1_3subset_variation_0.95<-bind_rows(sc1_3subset_variation_0.95)

 sc1_3subset_variation_0.95<-sc1_3subset_variation_0.95 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_3subset_variation_0.95)<-c("C",
                                         "q1","q2","q3",
                                         "pa1","pa2","pa3",
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
 
 sc1_3subset_variation_0.90<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,0.5,0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       1/3,1/3,1/3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,0.9,
       0.10,
       C_vars
       
     )
   }
   sc1_3subset_variation_0.90[[m]]<- data.frame(res)
 }
 saveRDS(sc1_3subset_variation_0.90,file="sc1_3subset_variation_estimate_0.90.rds")
 sc1_3subset_variation_0.90<-bind_rows(sc1_3subset_variation_0.90)
 
 sc1_3subset_variation_0.90<-sc1_3subset_variation_0.90 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_3subset_variation_0.90)<-c("C",
                                           "q1","q2","q3",
                                           "pa1","pa2","pa3",
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
 
 sc1_3subset_variation_0.85<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,0.5,0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       1/3,1/3,1/3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,0.9,
       0.15,
       C_vars
       
     )
   }
   sc1_3subset_variation_0.85[[m]]<- data.frame(res)
 }
 saveRDS(sc1_3subset_variation_0.85,file="sc1_3subset_variation_estimate_0.85.rds")
 sc1_3subset_variation_0.85<-bind_rows(sc1_3subset_variation_0.85)
 
 sc1_3subset_variation_0.85<-sc1_3subset_variation_0.85 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_3subset_variation_0.85)<-c("C",
                                           "q1","q2","q3",
                                           "pa1","pa2","pa3",
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
 
 
 sc1_3subset_variation_0.80<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,0.5,0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       0.4,
       1/3,1/3,1/3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,0.9,
       0.20,
       C_vars
       
     )
   }
   sc1_3subset_variation_0.80[[m]]<- data.frame(res)
 }
 saveRDS(sc1_3subset_variation_0.80,file="sc1_3subset_variation_estimate_0.80.rds")
 sc1_3subset_variation_0.80<-bind_rows(sc1_3subset_variation_0.80)
 
 sc1_3subset_variation_0.80<-sc1_3subset_variation_0.80 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_3subset_variation_0.80)<-c("C",
                                           "q1","q2","q3",
                                           "pa1","pa2","pa3",
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


 ####Scenarios with 2 subset and under the null hypothesis ####
 
 sc1_2subset_variation_0.80<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]

   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,0.5,
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,0.9,
       0.20,
       C_vars

     )
   }
   sc1_2subset_variation_0.80[[m]]<- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_0.80,file="sc1_2subset_variation_estimate_0.80.rds")
 sc1_2subset_variation_0.80<-bind_rows(sc1_2subset_variation_0.80)

 sc1_2subset_variation_0.80<-sc1_2subset_variation_0.80 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_0.80)<-c("C",
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
                                         "interaction_quali_1",
                                         "interaction_quanti_1",
                                         "enrichment0_1",
                                         "enrichment1_1",
                                         "continue_2 with subset ",
                                         "efficacy0_2",
                                         "efficacy1_2",
                                         "interaction_quali_2",
                                         "interaction_quanti_2",
                                         "enrichment0_2",
                                         "enrichment1_2",
                                         "continue_3 with subset ",
                                         "efficacy0_3",
                                         "efficacy1_3",
                                         "interaction_quali_3",
                                         "interaction_quanti_3",
                                         "enrichment0_3",
                                         "enrichment1_3",
                                         "continue_4 with subset ",
                                         "efficacy0_4",
                                         "efficacy1_4",
                                         "interaction_quali_4",
                                         "interaction_quanti_4",
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
                                         "biais_1_sd"
 )
 sc1_2subset_variation_0.85<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,0.5,
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,0.9,
       0.15,
       C_vars
       
     )
   }
   sc1_2subset_variation_0.80[[m]]<- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_0.85,file="sc1_2subset_variation_estimate_0.85.rds")
 sc1_2subset_variation_0.85<-bind_rows(sc1_2subset_variation_0.85)
 
 sc1_2subset_variation_0.85<-sc1_2subset_variation_0.85 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_0.85)<-c("C",
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
                                         "interaction_quali_1",
                                         "interaction_quanti_1",
                                         "enrichment0_1",
                                         "enrichment1_1",
                                         "continue_2 with subset ",
                                         "efficacy0_2",
                                         "efficacy1_2",
                                         "interaction_quali_2",
                                         "interaction_quanti_2",
                                         "enrichment0_2",
                                         "enrichment1_2",
                                         "continue_3 with subset ",
                                         "efficacy0_3",
                                         "efficacy1_3",
                                         "interaction_quali_3",
                                         "interaction_quanti_3",
                                         "enrichment0_3",
                                         "enrichment1_3",
                                         "continue_4 with subset ",
                                         "efficacy0_4",
                                         "efficacy1_4",
                                         "interaction_quali_4",
                                         "interaction_quanti_4",
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
                                         "biais_1_sd"
 )
 
 sc1_2subset_variation_0.90<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,0.5,
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,0.9,
       0.10,
       C_vars
       
     )
   }
   sc1_2subset_variation_0.90[[m]]<- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_0.90,file="sc1_2subset_variation_estimate_0.90.rds")
 sc1_2subset_variation_0.90<-bind_rows(sc1_2subset_variation_0.90)
 
 sc1_2subset_variation_0.90<-sc1_2subset_variation_0.90 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_0.90)<-c("C",
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
                                         "interaction_quali_1",
                                         "interaction_quanti_1",
                                         "enrichment0_1",
                                         "enrichment1_1",
                                         "continue_2 with subset ",
                                         "efficacy0_2",
                                         "efficacy1_2",
                                         "interaction_quali_2",
                                         "interaction_quanti_2",
                                         "enrichment0_2",
                                         "enrichment1_2",
                                         "continue_3 with subset ",
                                         "efficacy0_3",
                                         "efficacy1_3",
                                         "interaction_quali_3",
                                         "interaction_quanti_3",
                                         "enrichment0_3",
                                         "enrichment1_3",
                                         "continue_4 with subset ",
                                         "efficacy0_4",
                                         "efficacy1_4",
                                         "interaction_quali_4",
                                         "interaction_quanti_4",
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
                                         "biais_1_sd"
 )
 
 sc1_2subset_variation_0.95<-list()
 for(m in seq_along(C_var)){
   C_vars<-C_var[m]
   
   res<- foreach(i=1:Nsimu,
                 .inorder=F,.combine="rbind",
                 .packages=c("dplyr","R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,0.5,
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,0.9,
       0.05,
       C_vars
       
     )
   }
   sc1_2subset_variation_0.95[[m]]<- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_0.95,file="sc1_2subset_variation_estimate_0.95.rds")
 sc1_2subset_variation_0.95<-bind_rows(sc1_2subset_variation_0.95)
 
 sc1_2subset_variation_0.95<-sc1_2subset_variation_0.95 %>% group_by(C) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_0.95)<-c("C",
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
                                         "interaction_quali_1",
                                         "interaction_quanti_1",
                                         "enrichment0_1",
                                         "enrichment1_1",
                                         "continue_2 with subset ",
                                         "efficacy0_2",
                                         "efficacy1_2",
                                         "interaction_quali_2",
                                         "interaction_quanti_2",
                                         "enrichment0_2",
                                         "enrichment1_2",
                                         "continue_3 with subset ",
                                         "efficacy0_3",
                                         "efficacy1_3",
                                         "interaction_quali_3",
                                         "interaction_quanti_3",
                                         "enrichment0_3",
                                         "enrichment1_3",
                                         "continue_4 with subset ",
                                         "efficacy0_4",
                                         "efficacy1_4",
                                         "interaction_quali_4",
                                         "interaction_quanti_4",
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
                                         "biais_1_sd"
 )

 saveRDS(sc1_3subset_variation_0.95, file = "sc1_3subset_variation_final_0.95.rds")
 saveRDS(sc1_3subset_variation_0.90, file = "sc1_3subset_variation_final_0.90.rds")
 saveRDS(sc1_3subset_variation_0.85, file = "sc1_3subset_variation_final_0.85.rds")
 saveRDS(sc1_3subset_variation_0.80, file = "sc1_3subset_variation_final_0.80.rds")
 
 saveRDS(sc1_2subset_variation_0.95, file = "sc1_2subset_variation_final_0.95.rds")
 saveRDS(sc1_2subset_variation_0.90, file = "sc1_2subset_variation_final_0.90.rds")
 saveRDS(sc1_2subset_variation_0.85, file = "sc1_2subset_variation_final_0.85.rds")
 saveRDS(sc1_2subset_variation_0.80, file = "sc1_2subset_variation_final_0.80.rds")


 #### data management for graphs #####
 sc1_3_subset_variation<-rbind(sc1_3subset_variation_final_0.80,sc1_3subset_variation_final_0.85,sc1_3subset_variation_final_0.90,sc1_3subset_variation_final_0.95)
 threshold<-rep(c(0.20,0.15,0.10,0.05),each=401)      
 threshold<-c(threshold)
 sc1_3_subset_variation<-cbind(sc1_3_subset_variation,threshold)
 sc1_3_subset_variation$threshold<-as.factor(sc1_3_subset_variation$threshold)
 
 sc1_2_subset_variation<-rbind(sc1_2subset_variation_final_0.80,sc1_2subset_variation_final_0.85,sc1_2subset_variation_final_0.90,sc1_2subset_variation_final_0.95)
 threshold<-rep(c(0.20,0.15,0.10,0.05),each=401)      
 threshold<-c(threshold)
 sc1_2_subset_variation<-cbind(sc1_2_subset_variation,threshold)
 sc1_2_subset_variation$threshold<-as.factor(sc1_2_subset_variation$threshold)
 
 ###graphs 
 ggplot(
   data.frame(sc1_3_subset_variation),
   aes(
     x = sc1_3_subset_variation$interaction_quali_4,
     y = C,
     fill = threshold
   )
 ) + geom_line(aes(color = threshold)) + annotate(
   geom = "point",
   x = c(0.20, 0.15, 0.10, 0.05),
   y = c(0.325, 0.505, 0.980, 1.970),
   shape = 4,
   size = 4
 ) + scale_y_continuous(breaks =   seq(0, 2, 0.1)) + scale_x_reverse(breaks =   seq(0, 1, 0.05)) +
   xlab("Interaction qualitative proportion") + ylab("Critical values") + labs(color =
                                                                                 "Threshold") + theme_classic() + theme(legend.background = element_rect(
                                                                                   size=0.5, linetype="solid", 
                                                                                   colour ="black"),
                                                                                   legend.title = element_text(size = 18),
                                                                                   legend.text = element_text(size = 10),
                                                                                   axis.text.x = element_text(size = 10),
                                                                                   axis.text.y = element_text(size = 10) ,
                                                                                   axis.title.y = element_text(size = 18),
                                                                                   axis.title.x = element_text(size = 18)
                                                                                 ) +geom_hline(yintercept = c(0.325, 0.505, 0.980, 1.970),linetype="dashed",alpha=0.4)+scale_color_npg()
 ggplot(
   data.frame(sc1_3_subset_variation),
   aes(
     x = sc1_3_subset_variation$interaction_quanti_4,
     y = C,
     fill = threshold
   )
 ) + geom_line(aes(color = threshold)) + annotate(
   geom = "point",
   x = c(0.20, 0.15, 0.10, 0.05),
   y = c(0.155, 0.230, 0.360, 0.710),
   shape = 4,
   size = 4
 ) + scale_y_continuous(breaks =   seq(0, 1, 0.1)) + scale_x_reverse(breaks =   seq(0, 1, 0.05)) +
   xlab("Interaction quantitative proportion") + ylab("Critical values") + labs(color =
                                                                                  "Threshold") + theme_classic() + theme(legend.background = element_rect(
                                                                                    size=0.5, linetype="solid", 
                                                                                    colour ="black"),
                                                                                    legend.title = element_text(size = 18),
                                                                                    legend.text = element_text(size = 10),
                                                                                    axis.text.x = element_text(size = 10),
                                                                                    axis.text.y = element_text(size = 10) ,
                                                                                    axis.title.y = element_text(size = 18),
                                                                                    axis.title.x = element_text(size = 18)
                                                                                  ) +geom_hline(yintercept = c(0.155, 0.230, 0.360, 0.710),linetype="dashed",alpha=0.4)+scale_color_npg()

 
 
 
 
 ggplot(
   data.frame(sc1_2_subset_variation),
   aes(
     x = sc1_2_subset_variation$interaction_quanti_4,
     y = C,
     fill = threshold
   )
 ) + geom_line(aes(color = threshold)) + annotate(
   geom = "point",
   x = c(0.20, 0.15, 0.10, 0.05),
   y = c(0.045,0.070,0.105,0.220),
   shape = 4,
   size = 4
 ) + scale_y_continuous(breaks =   seq(0, 1, 0.1)) + scale_x_reverse(breaks =   seq(0, 1, 0.05)) +
   xlab("Interaction quantitative proportion") + ylab("Critical values") + labs(color =
                                                                                  "Threshold") + theme_classic() + theme(legend.background = element_rect(
                                                                                    size=0.5, linetype="solid", 
                                                                                    colour ="black"),
                                                                                    legend.title = element_text(size = 18),
                                                                                    legend.text = element_text(size = 10),
                                                                                    axis.text.x = element_text(size = 10),
                                                                                    axis.text.y = element_text(size = 10) ,
                                                                                    axis.title.y = element_text(size = 18),
                                                                                    axis.title.x = element_text(size = 18)
                                                                                  ) +geom_hline(yintercept = c(0.045,0.070,0.105,0.220),linetype="dashed",alpha=0.4)+scale_color_npg()
 
 ggplot(
   data.frame(sc1_2_subset_variation),
   aes(
     x = sc1_2_subset_variation$interaction_quali_4,
     y = C,
     fill = threshold
   )
 ) + geom_line(aes(color = threshold)) + annotate(
   geom = "point",
   x = c(0.20, 0.15, 0.10, 0.05),
   y = c(0.060,0.110,0.210,0.655),
   shape = 4,
   size = 4
 ) + scale_y_continuous(breaks =   seq(0, 1, 0.1)) + scale_x_reverse(breaks =   seq(0, 1, 0.05)) +
   xlab("Interaction qualitative proportion") + ylab("Critical values") + labs(color =
                                                                                 "Threshold") + theme_classic() + theme(legend.background = element_rect(
                                                                                   size=0.5, linetype="solid", 
                                                                                   colour ="black"),
                                                                                   legend.title = element_text(size = 18),
                                                                                   legend.text = element_text(size = 10),
                                                                                   axis.text.x = element_text(size = 10),
                                                                                   axis.text.y = element_text(size = 10) ,
                                                                                   axis.title.y = element_text(size = 18),
                                                                                   axis.title.x = element_text(size = 18)
                                                                                 ) +geom_hline(yintercept =  c(0.060,0.110,0.210,0.655),linetype="dashed",alpha=0.4)+scale_color_npg()
 