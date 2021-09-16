#Load packages
install.packages("R2jags")
install.packages("dplyr")
install.packages("doSNOW")
library(R2jags)
library(dplyr)
library(parallel)
library(foreach)
library(doSNOW)

### values explication 
#for 2 subsets 
#pi1A= death proportion in subset A and group 1
#pi1B= death proportion in subset B and group 1
#pi0A= death proportion in subset A and group 0
#pi0B= death proportion in subset B and group 0
# for 3 subsets
#pi1A= death proportion in subset A and group 1
#pi1A= death proportion in subset A and group 0
#pi1B= death proportion in subset B and group 1
#pi0B= death proportion in subset B and group 0
#pi1C= death proportion in subset C and group 1
#pi0C= death proportion in subset C and group 0

#lambda_1 = clinical relevant threshold for efficacy measure
# C1 and C2 are critical values for GS interaction test ( C1 qualitative and C2 quantitative )
#gamma_1 = decision threshold for efficacy measure
#espilon_1 = decision threshold fo interaction measure 
#Q1 = balance of randomization
#pa = prevalence of subset A
#RR.subset 1 = relative risk in subset A
#RR.subset 0 = relative risk in subset B
#N_1,N_2,N_3,N_4 = sample size recruitment at each interim and terminal analysis

### load function for 3susbets######
####subset 3####

### stopping rules function for first analysis#####
stopping_rules_first_analyse_variation_3_subset <-
  function(P1_1, P1_2, P1_3, P5, P6, gamma, epsilon) {
    nDecision <- "continue with 3 subsets"
    continue <- 0
    enrichment_1 <- 0
    enrichment_2 <- 0
    enrichment_3 <- 0
    efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
    efficacy_2 <- ifelse(P1_2 > gamma, 1, 0)
    efficacy_3 <- ifelse(P1_3 > gamma, 1, 0)
    interaction_quali <- ifelse(P5 > epsilon, 1, 0)
    interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
    interaction_type <-
      ifelse(
        interaction_quali == 1 & interaction_quanti == 1,
        "both interaction",
        ifelse(
          interaction_quali == 1,
          "quali",
          ifelse(interaction_quanti == 1, "quanti", "NO")
        )
      )
    enrichment_1 <-
      ifelse(efficacy_1 == 1 &
               (interaction_quali == 1 | interaction_quanti == 1),
             1,
             0)
    enrichment_2 <-
      ifelse(efficacy_2 == 1 &
               (interaction_quali == 1 | interaction_quanti  == 1),
             1,
             0)
    enrichment_3 <-
      ifelse(efficacy_3 == 1 &
               (interaction_quali == 1 | interaction_quanti  == 1),
             1,
             0)
    enrichment_1_3 <-
      ifelse(efficacy_1 == 1 &
               efficacy_3 == 1 &
               (interaction_quali == 1 | interaction_quanti  == 1),
             1,
             0)
    enrichment_1_2 <-
      ifelse(efficacy_1 == 1 &
               efficacy_2 == 1 &
               (interaction_quali == 1 | interaction_quanti  == 1),
             1,
             0)
    enrichment_2_3 <-
      ifelse(efficacy_2 == 1 &
               efficacy_3 == 1 &
               (interaction_quali == 1 | interaction_quanti  == 1),
             1,
             0)
    
    if (enrichment_1 == 1 &
        enrichment_2 == 1 & enrichment_3 == 0) {
      nDecision <- "efficacy in the subset 1 & 2"
      enrichment_1 <- 0
      efficacy_1 <- 1
      enrichment_2 <- 0
      efficacy_2 <- 1
      enrichment_1_2 <- 1
      continue <- 0
    }
    else if (enrichment_1 == 1 &
             enrichment_2 == 0 & enrichment_3 == 1) {
      nDecision <- "efficacy in the subset 1 & 3"
      enrichment_1 <- 0
      efficacy_1 <- 1
      efficacy_3 <- 1
      enrichment_3 <- 0
      enrichment_1_3 <- 1
      continue <- 0
    }
    else if (enrichment_1 == 0 &
             enrichment_2 == 1 & enrichment_3 == 1) {
      nDecision <- "efficacy in the subset 2 & 3"
      enrichment_2 <- 0
      efficacy_2 <- 1
      efficacy_3 <- 1
      enrichment_3 <- 0
      enrichment_2_3 <- 1
      continue <- 0
    }
    else if (enrichment_1 == 1 &
             enrichment_2 == 0 & enrichment_3 == 0) {
      nDecision <- "efficacy in the subset 1"
      continue <- 0
      enrichment_1 <- 1
      efficacy_1 <- 1
      
    }
    else if (enrichment_1 == 0 &
             enrichment_2 == 1 & enrichment_3 == 0) {
      nDecision <- "efficacy in the subset 2"
      continue <- 0
      enrichment_2 <- 1
      efficacy_2 <- 1
      
    }
    else if (enrichment_1 == 0 &
             enrichment_2 == 0 & enrichment_3 == 1) {
      nDecision <- "efficacy in the subset 3"
      continue <- 0
      enrichment_3 <- 1
      efficacy_3 <- 1
      
    }
    else{
      nDecision <- "continue with 3 subsets"
      continue <- 1
      enrichment_1 <- 0
      enrichment_2 <- 0
      enrichment_3 <- 0
      enrichment_1_2 <- 0
      enrichment_1_3 <- 0
      enrichment_2_3 <- 0
    }
    return(
      c(
        nDecision,
        continue,
        efficacy_1,
        efficacy_2,
        efficacy_3,
        interaction_quali,
        interaction_quanti,
        enrichment_1,
        enrichment_2,
        enrichment_3,
        enrichment_1_2,
        enrichment_1_3,
        enrichment_2_3,
        interaction_type
      )
    )
  }


### stopping rules function for 2,3 and terminal analysis#####
stopping_rules_other_analyse_variation_3_subset <-
  function(nDecision,
           P1_1,
           P1_2,
           P1_3,
           P5,
           P6,
           gamma,
           epsilon,
           interaction_type) {
    continue <- 0
    efficacy_1 <- 0
    efficacy_2 <- 0
    efficacy_3 <- 0
    interaction_quali <- 0
    interaction_quanti <- 0
    enrichment_1 <- 0
    enrichment_2 <- 0
    enrichment_3 <- 0
    enrichment_1_2 <- 0
    enrichment_1_3 <- 0
    enrichment_2_3 <- 0
    
    if (nDecision == "efficacy in the subset 1 & 2") {
      efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
      efficacy_2 <- ifelse(P1_2 > gamma, 1, 0)
      efficacy_3 <- ifelse(P1_3 > gamma, 1, 0)
      interaction_quali <- ifelse(P5 > epsilon, 1, 0)
      interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
      interaction_type <-
        ifelse(
          interaction_quali == 1 & interaction_quanti == 1,
          "both interaction",
          ifelse(
            interaction_quali == 1,
            "quali",
            ifelse(interaction_quanti == 1, "quanti", "NO")
          )
        )
      enrichment_1 <-
        ifelse(efficacy_1 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_2 <-
        ifelse(efficacy_2 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_3 <-
        ifelse(efficacy_3 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_1_2 <-
        ifelse(
          efficacy_1 == 1 &
            efficacy_2 == 1 &
            (interaction_quali == 1 |
               interaction_quanti  == 1),
          1,
          0
        )
      
      if (enrichment_1 == 1 &
          enrichment_2 == 1 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 1 & 2"
        enrichment_1 <- 0
        efficacy_1 <- 1
        enrichment_2 <- 0
        efficacy_2 <- 1
        enrichment_1_2 <- 1
        continue <- 0
      }
      else if (enrichment_1 == 1 &
               enrichment_2 == 0 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 1"
        continue <- 0
        enrichment_1 <- 1
        efficacy_1 <- 1
        
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 1 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 2"
        continue <- 0
        enrichment_2 <- 1
        efficacy_2 <- 1
        
      }
    } else if (nDecision == "efficacy in the subset 1 & 3") {
      efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
      efficacy_2 <- ifelse(P1_2 > gamma, 1, 0)
      efficacy_3 <- ifelse(P1_3 > gamma, 1, 0)
      interaction_quali <- ifelse(P5 > epsilon, 1, 0)
      interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
      interaction_type <-
        ifelse(
          interaction_quali == 1 & interaction_quanti == 1,
          "both interaction",
          ifelse(
            interaction_quali == 1,
            "quali",
            ifelse(interaction_quanti == 1, "quanti", "NO")
          )
        )
      enrichment_1 <-
        ifelse(efficacy_1 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_2 <-
        ifelse(efficacy_2 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_3 <-
        ifelse(efficacy_3 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_1_3 <-
        ifelse(
          efficacy_1 == 1 &
            efficacy_3 == 1 &
            (interaction_quali == 1 |
               interaction_quanti  == 1),
          1,
          0
        )
      
      if (enrichment_1 == 1 &
          enrichment_2 == 0 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 1 & 3"
        enrichment_1 <- 0
        efficacy_1 <- 1
        efficacy_3 <- 1
        enrichment_3 <- 0
        enrichment_1_3 <- 1
        continue <- 0
      }
      
      else if (enrichment_1 == 1 &
               enrichment_2 == 0 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 1"
        continue <- 0
        enrichment_1 <- 1
        efficacy_1 <- 1
        
      }
      
      else if (enrichment_1 == 0 &
               enrichment_2 == 0 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 3"
        continue <- 0
        enrichment_3 <- 1
        efficacy_3 <- 1
        
      }
    } else if (nDecision == "efficacy in the subset 2 & 3") {
      efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
      efficacy_2 <- ifelse(P1_2 > gamma, 1, 0)
      efficacy_3 <- ifelse(P1_3 > gamma, 1, 0)
      interaction_quali <- ifelse(P5 > epsilon, 1, 0)
      interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
      interaction_type <-
        ifelse(
          interaction_quali == 1 & interaction_quanti == 1,
          "both interaction",
          ifelse(
            interaction_quali == 1,
            "quali",
            ifelse(interaction_quanti == 1, "quanti", "NO")
          )
        )
      enrichment_1 <-
        ifelse(efficacy_1 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_2 <-
        ifelse(efficacy_2 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_3 <-
        ifelse(efficacy_3 == 1 &
                 (interaction_quali == 1 |
                    interaction_quanti  == 1),
               1,
               0)
      enrichment_2_3 <-
        ifelse(
          efficacy_2 == 1 &
            efficacy_3 == 1 &
            (interaction_quali == 1 |
               interaction_quanti  == 1),
          1,
          0
        )
      
      if (enrichment_1 == 0 &
          enrichment_2 == 1 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 2 & 3"
        enrichment_2 <- 0
        efficacy_2 <- 1
        efficacy_3 <- 1
        enrichment_3 <- 0
        enrichment_2_3 <- 1
        continue <- 0
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 1 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 2"
        continue <- 0
        enrichment_2 <- 1
        efficacy_2 <- 1
        
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 0 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 3"
        continue <- 0
        enrichment_3 <- 1
        efficacy_3 <- 1
        
      }
    } else if (nDecision == "efficacy in the subset 1") {
      nDecision <- "efficacy in the subset 1"
      enrichment_1 <- 1
      efficacy_1 <- 1
      interaction_quali <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quali"
        ),
        1,
        0)
      interaction_quanti <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quanti"
        ),
        1,
        0)
      continue <- 0
    }
    else if (nDecision == "efficacy in the subset 2") {
      nDecision <- "efficacy in the subset 2"
      enrichment_2 <- 1
      efficacy_2 <- 1
      interaction_quali <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quali"
        ),
        1,
        0)
      interaction_quanti <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quanti"
        ),
        1,
        0)
      continue <- 0
    }
    else if (nDecision == "efficacy in the subset 3") {
      nDecision <- "efficacy in the subset 3"
      enrichment_3 <- 1
      efficacy_3 <- 1
      interaction_quali <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quali"
        ),
        1,
        0)
      interaction_quanti <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quanti"
        ),
        1,
        0)
      continue <- 0
    }
    else {
      continue <- 0
      efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
      efficacy_2 <- ifelse(P1_2 > gamma, 1, 0)
      efficacy_3 <- ifelse(P1_3 > gamma, 1, 0)
      interaction_quali <- ifelse(P5 > epsilon, 1, 0)
      interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
      interaction_type <-
        ifelse(
          interaction_quali == 1 & interaction_quanti == 1,
          "both interaction",
          ifelse(
            interaction_quali == 1,
            "quali",
            ifelse(interaction_quanti == 1, "quanti", "NO")
          )
        )
      enrichment_1 <-
        ifelse(efficacy_1 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_2 <-
        ifelse(efficacy_2 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_3 <-
        ifelse(efficacy_3 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_1_3 <-
        ifelse(
          efficacy_1 == 1 &
            efficacy_3 == 1 &
            (interaction_quali == 1 | interaction_quanti == 1),
          1,
          0
        )
      enrichment_1_2 <-
        ifelse(
          efficacy_1 == 1 &
            efficacy_2 == 1 &
            (interaction_quali == 1 | interaction_quanti == 1),
          1,
          0
        )
      enrichment_2_3 <-
        ifelse(
          efficacy_2 == 1 &
            efficacy_3 == 1 &
            (interaction_quali == 1 | interaction_quanti == 1),
          1,
          0
        )
      
      if (enrichment_1 == 1 &
          enrichment_2 == 1 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 1 & 2"
        enrichment_1 <- 0
        efficacy_1 <- 1
        enrichment_2 <- 0
        efficacy_2 <- 1
        enrichment_1_2 <- 1
        continue <- 0
      }
      else if (enrichment_1 == 1 &
               enrichment_2 == 0 & enrichment_3 == 3) {
        nDecision <- "efficacy in the subset 1 & 3"
        enrichment_1 <- 0
        efficacy_1 <- 1
        enrichment_3 <- 0
        efficacy_3 <- 1
        enrichment_1_3 <- 1
        continue <- 0
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 1 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 2 & 3"
        enrichment_2 <- 0
        efficacy_2 <- 1
        enrichment_3 <- 0
        efficacy_3 <- 1
        enrichment_2_3 <- 1
        continue <- 0
      }
      else if (enrichment_1 == 1 &
               enrichment_2 == 0 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 1"
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 1 & enrichment_3 == 0) {
        nDecision <- "efficacy in the subset 2"
      }
      else if (enrichment_1 == 0 &
               enrichment_2 == 0 & enrichment_3 == 1) {
        nDecision <- "efficacy in the subset 3"
      }
      else{
        nDecision <- "continue with 3 subsets"
        continue <- 1
        enrichment_1 <- 0
        enrichment_2 <- 0
        enrichment_3 <- 0
        enrichment_1_2 <- 0
        enrichment_1_3 <- 0
        enrichment_2_3 <- 0
      }
    }
    return(
      c(
        nDecision,
        continue,
        efficacy_1,
        efficacy_2,
        efficacy_3,
        interaction_quali,
        interaction_quanti,
        enrichment_1,
        enrichment_2,
        enrichment_3,
        enrichment_1_2,
        enrichment_1_3,
        enrichment_2_3,
        interaction_type
      )
    )
  }

### Matrix function for method with 3 subsets ####
create_matrix_3_subset <- function(nrow) {
  mat.2 <-
    matrix(nrow = nrow,
           ncol = 87,
           dimnames = list(
             NULL,
             c(
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
           ))
  
} #create matrix for simulation

##### functions for data simulation with 3 susbets ( each case )
create_data_frame_for_continue_with_3_subsets <- function(q1,
                                                          q2,
                                                          q3,
                                                          N_patients,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B,
                                                          pi1C,
                                                          pi0C,
                                                          pa1,
                                                          pa2,
                                                          pa3) {
  subset <-
    sample(1:3,
           N_patients,
           replace = T,
           prob = c(pa1, pa2, pa3))  #subsets
   #randomization
  group <- NULL
  group[subset == 1] <- rbinom(sum(subset == 1), 1, q1)
  group[subset == 2] <- rbinom(sum(subset == 2), 1, q2)
  group[subset == 3] <-
    rbinom(sum(subset == 3), 1, q3)  #Group traitement vs Placebo
   #Subset biomarkers
  y <- NULL                  #Response
  y[group == 1 & subset == 1] <-
    rbinom(sum(group == 1 & subset == 1), 1, pi1A)
  y[group == 0 & subset == 1] <-
    rbinom(sum(group == 0 & subset == 1), 1, pi0A)
  y[group == 1 & subset == 2] <-
    rbinom(sum(group == 1 & subset == 2), 1, pi1B)
  y[group == 0 & subset == 2] <-
    rbinom(sum(group == 0 & subset == 2), 1, pi0B)
  y[group == 1 & subset == 3] <-
    rbinom(sum(group == 1 & subset == 3), 1, pi1C)
  y[group == 0 & subset == 3] <-
    rbinom(sum(group == 0 & subset == 3), 1, pi0C)
  
  return(data.frame(group, subset, y))
}  #create data if decision rules aren't met


create_data_frame_if_efficacy_subset1_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <- NULL
    subset <- c(rep(1, N_patients))
    group <- NULL
    group[subset == 1] <-
      rbinom(sum(subset == 1), 1, q1)  #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                 # Response
    y[group == 1 &
        subset == 1] <-
      rbinom(sum(group == 1 & subset == 1), 1, pi1A)
    y[group == 0 &
        subset == 1] <-
      rbinom(sum(group == 0 & subset == 1), 1, pi0A)
    return(data.frame(group, subset, y))
  }  #create data if subset 1 is effective

create_data_frame_if_efficacy_subset2_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <- NULL
    subset <- c(rep(2, N_patients))
    group <- NULL
    group[subset == 2] <-
      rbinom(sum(subset == 2), 1, q2)  #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                  #Response
    y[group == 1 &
        subset == 2] <-
      rbinom(sum(group == 1 & subset == 2), 1, pi1B)
    y[group == 0 &
        subset == 2] <-
      rbinom(sum(group == 0 & subset == 2), 1, pi0B)
    return(data.frame(group, subset, y))
  }  #create data if subset 2 is effective

create_data_frame_if_efficacy_subset3_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <- NULL
    subset <- c(rep(3, N_patients))
    group <- NULL
    group[subset == 3] <-
      rbinom(sum(subset == 3), 1, q3)  #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                 # Response
    y[group == 1 &
        subset == 3] <-
      rbinom(sum(group == 1 & subset == 3), 1, pi1C)
    y[group == 0 &
        subset == 3] <-
      rbinom(sum(group == 0 & subset == 3), 1, pi0C)
    return(data.frame(group, subset, y))
  }  #create data if subset 3 is effective

create_data_frame_if_efficacy_subset1.2_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <-
      sample(c(1, 2),
             N_patients,
             replace = T,
             prob = c(pa1, pa2))  #subsets
     #randomization
    group <- NULL
    group[subset == 1] <- rbinom(sum(subset == 1), 1, q1)
    group[subset == 2] <-
      rbinom(sum(subset == 2), 1, q2) #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                  #Response
    y[group == 1 & subset == 1] <-
      rbinom(sum(group == 1 & subset == 1), 1, pi1A)
    y[group == 0 & subset == 1] <-
      rbinom(sum(group == 0 & subset == 1), 1, pi0A)
    y[group == 1 & subset == 2] <-
      rbinom(sum(group == 1 & subset == 2), 1, pi1B)
    y[group == 0 & subset == 2] <-
      rbinom(sum(group == 0 & subset == 2), 1, pi0B)
    
    return(data.frame(group, subset, y))
  }  #create data if subset 1/2 is effective

create_data_frame_if_efficacy_subset1.3_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <-
      sample(c(1, 3),
             N_patients,
             replace = T,
             prob = c(pa1, pa3)) # subsets
     #randomization
    group <- NULL
    group[subset == 1] <- rbinom(sum(subset == 1), 1, q1)
    group[subset == 3] <-
      rbinom(sum(subset == 3), 1, q3) #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                 # Response
    y[group == 1 & subset == 1] <-
      rbinom(sum(group == 1 & subset == 1), 1, pi1A)
    y[group == 0 & subset == 1] <-
      rbinom(sum(group == 0 & subset == 1), 1, pi0A)
    y[group == 1 & subset == 3] <-
      rbinom(sum(group == 1 & subset == 3), 1, pi1C)
    y[group == 0 & subset == 3] <-
      rbinom(sum(group == 0 & subset == 3), 1, pi0C)
    
    return(data.frame(group, subset, y))
  }  #create data if subset 1/3 is effective

create_data_frame_if_efficacy_subset2.3_3_subsets <-
  function(q1,
           q2,
           q3,
           N_patients,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3) {
    subset <-
      sample(c(2, 3),
             N_patients,
             replace = T,
             prob = c(pa2, pa3)) # subsets
     #randomization
    group <- NULL
    group[subset == 2] <- rbinom(sum(subset == 2), 1, q2)
    group[subset == 3] <-
      rbinom(sum(subset == 3), 1, q3) #Group traitement vs Placebo
     #Subset biomarkers
    y <- NULL                 # Response
    y[group == 1 & subset == 2] <-
      rbinom(sum(group == 1 & subset == 2), 1, pi1B)
    y[group == 0 & subset == 2] <-
      rbinom(sum(group == 0 & subset == 2), 1, pi0B)
    y[group == 1 & subset == 3] <-
      rbinom(sum(group == 1 & subset == 3), 1, pi1C)
    y[group == 0 & subset == 3] <-
      rbinom(sum(group == 0 & subset == 3), 1, pi0C)
    
    return(data.frame(group, subset, y))
  }  #create data if subset 2/3 is effective


#### Bayesian model for JAGS and R2jags
binary_bayesian_model_3_subset <- function() {
  nb11 ~ dbin(p11, n11)
  nb12 ~ dbin(p12, n12)
  nb13 ~ dbin(p13, n13)
  nb01 ~ dbin(p01, n01)
  nb02 ~ dbin(p02, n02)
  nb03 ~ dbin(p03, n03)
  
  p11 ~ dbeta(1, 1)
  p01 ~ dbeta(1, 1)
  p12 ~ dbeta(1, 1)
  p02 ~ dbeta(1, 1)
  p13 ~ dbeta(1, 1)
  p03 ~ dbeta(1, 1)
  
  
  p1.estim <-
    ((p11 * pa1) + (p12 * pa2) + (p13 * pa3)) # probability in treatment group
  p0.estim <-
    ((p01 * pa1) + (p02 * pa2) + (p03 * pa3))  # probability in control group
  RR <- p1.estim / p0.estim
  RR.subset1 <- p11 / p01  #effect of treatment in subset 1
  RR.subset2 <- p12 / p02  #effect of treatment in subset 2
  RR.subset3 <- p13 / p03  #effect of treatment in subset 3
  beta_1 <- log(RR.subset1)
  beta_2 <- log(RR.subset2)
  beta_3 <- log(RR.subset3)
  sum1 <- (RR.subset1 - RR)
  sum2 <- (RR.subset2 - RR)
  sum3 <- (RR.subset3 - RR)
  sum_final <- sum(sum1, sum2, sum3)
  P6 <- (sum_final > C2)
  sd1 <- sqrt(mean((y1 - mean(y1)) ^ 2))
  sd2 <- sqrt(mean((y2 - mean(y2)) ^ 2))
  sd3 <- sqrt(mean((y3 - mean(y3)) ^ 2))
  betas <- c(beta_1, beta_2, beta_3)
  sigmas <- c(sd1, sd2, sd3)
  over <- betas / sigmas
  oversq <- over ^ 2
  test_quali <-
    min(sum(oversq * (over > 0)), sum(oversq * (over < 0)))
  P5 <- (test_quali > C1)
  
  P1_1 <- (RR.subset1 < lambda)
  P1_2 <- (RR.subset2 < lambda)
  P1_3 <- (RR.subset3 < lambda)
  
  biais_1 <- RR.subset1 - (pi1A / pi0A)
  biais_2 <- RR.subset2 - (pi1B / pi0B)
  biais_3 <- RR.subset3 - (pi1C / pi0C)
  
} # bayesian model with 3 subset for jags function

#### function for initialization values for R2jags
create_list_for_inits_values_for_3_subsets <-
  function(Data,
           pi1A,
           pi0A,
           pi1B,
           pi0B,
           pi1C,
           pi0C,
           pa1,
           pa2,
           pa3,
           lambda,
           C1,
           C2) {
    list(
      n11 = sum(Data$group == 1 & Data$subset == 1),
      n01 = sum(Data$group == 0 &
                  Data$subset == 1),
      n12 = sum(Data$group == 1 &
                  Data$subset == 2),
      n02 = sum(Data$group == 0 &
                  Data$subset == 2),
      n13 = sum(Data$group == 1 &
                  Data$subset == 3),
      n03 = sum(Data$group == 0 &
                  Data$subset == 3),
      nb11 = sum(Data$group == 1 & Data$subset == 1 & Data$y == 1,
                 na.rm = T),
      nb01 = sum(Data$group == 0 & Data$subset == 1 & Data$y == 1,
                 na.rm = T),
      nb12 = sum(Data$group == 1 & Data$subset == 2 & Data$y == 1,
                 na.rm = T),
      nb02 = sum(Data$group == 0 & Data$subset == 2 & Data$y == 1,
                 na.rm = T),
      nb13 = sum(Data$group == 1 & Data$subset == 3 & Data$y == 1,
                 na.rm = T),
      nb03 = sum(Data$group == 0 & Data$subset == 3 & Data$y == 1,
                 na.rm = T),
      y1 = filter(Data, subset == 1)[, 3],
      y2 = filter(Data, subset == 2)[, 3],
      y3 = filter(Data, subset == 3)[, 3],
      pi1A = pi1A,
      pi0A = pi0A,
      pi1B = pi1B,
      pi0B = pi0B,
      pi1C = pi1C,
      pi0C = pi0C,
      pa1 = pa1,
      pa2 = pa2,
      pa3 = pa3,
      lambda = lambda,
      C1 = C1,
      C2 = C2
    )
  } #list for iniatialization values for jags


#### function for R2JAGS result 
results_bayesian_model_3_subsets <-
  function(bayesian_jags_fonction_IA1, Data) {
    p11 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p11
    p10 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p10
    p01 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p01
    p00 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$p00
    RR <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR
    RR_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[6, 3]
    RR_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[6, 7]
    RR.subset1 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR.subset1
    RR.subset1_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[7, 3]
    RR.subset1_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[7, 7]
    RR.subset2 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR.subset2
    RR.subset2_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[8, 3]
    RR.subset2_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[8, 7]
    RR.subset3 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$RR.subset3
    RR.subset3_2.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[9, 3]
    RR.subset3_97.5 <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[9, 7]
    beta_1 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$beta_1
    beta_1_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$sd$beta_1
    beta_2 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$beta_2
    beta_2_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$sd$beta_2
    beta_3 <-
      bayesian_jags_fonction_IA1$BUGSoutput$mean$beta_3
    beta_3_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$sd$beta_3
    P1_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P1_1
    P1_2 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P1_2
    P1_3 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P1_3
    P5 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P5
    P6 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P6
    
    
    biais_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_1
    biais_1_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[13, 2]
    biais_2 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_2
    biais_2_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[14, 2]
    biais_3 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_3
    biais_3_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[15, 2]
    
    N_subset1_stage1 <- sum(Data$subset == 1)
    N_subset2_stage1 <- sum(Data$subset == 2)
    N_subset3_stage1 <- sum(Data$subset == 3)
    results <-
      cbind(
        RR,
        RR_2.5,
        RR_97.5,
        RR.subset1,
        RR.subset1_2.5,
        RR.subset1_97.5,
        RR.subset2,
        RR.subset2_2.5,
        RR.subset2_97.5,
        RR.subset3,
        RR.subset3_2.5,
        RR.subset3_97.5,
        P1_1,
        P1_2,
        P1_3,
        P5,
        P6,
        beta_1,
        beta_1_sd,
        beta_2,
        beta_2_sd,
        beta_3,
        beta_3_sd,
        biais_1,
        biais_1_sd,
        biais_2,
        biais_2_sd,
        biais_3,
        biais_3_sd,
        N_subset1_stage1,
        N_subset2_stage1,
        N_subset3_stage1
      )
    results
  } #results from bayesian model


#### Simulation function for 3 subset method#####

simulate_binary_bayesian_trial_3_subsets <-
  function(q1 = 0.5,
           q2 = 0.5,
           q3 = 0.5,
           N_1 = 200,
           N_2 = 200,
           N_3 = 200,
           N_4 = 200,
           pi1A = 0.3,
           pi0A = 0.3,
           pi1B = 0.3,
           pi0B = 0.3,
           pi1C = 0.3,
           pi0C = 0.3,
           pa1 = 1 / 3,
           pa2 = 1 / 3,
           pa3 = 1 / 3,
           bayesian.model = binary_bayesian_model_3_subset,
           stopping_rules.IA = stopping_rules_first_analyse_variation_3_subset,
           stopping_rules.TA = stopping_rules_other_analyse_variation_3_subset,
           lambda = 0.9,
           gamma = 0.9,
           epsilon = 0.05,
           C1,
           C2) {
    mat.sub1 <- create_matrix_3_subset(Nsimu)
    #for (i in 1:Nsimu) {
    First_result <- NULL
    Second_result <- NULL
    Third_result <- NULL
    Final_result <- NULL
    
    Data_IA1 <-
      create_data_frame_for_continue_with_3_subsets(q1, q2, q3,
                                                    N_1,
                                                    pi1A,
                                                    pi0A,
                                                    pi1B,
                                                    pi0B, pi1C, pi0C,
                                                    pa1, pa2, pa3)
    
    inits_IA1 <-
      create_list_for_inits_values_for_3_subsets(Data_IA1,
                                                 pi1A,
                                                 pi0A,
                                                 pi1B,
                                                 pi0B,
                                                 pi1C,
                                                 pi0C,
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
    print(Check_interim_decision1)
    
    if (Check_interim_decision1[1] == "continue with 3 subsets") {
      Data_IA2 <- create_data_frame_for_continue_with_3_subsets(q1, q2, q3,
                                                                N_2,
                                                                pi1A,
                                                                pi0A,
                                                                pi1B,
                                                                pi0B, pi1C, pi0C,
                                                                pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset1.2_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset1.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset2.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset1_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset2_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA2 <-
        create_data_frame_if_efficacy_subset3_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      inits_IA2 <-
        create_list_for_inits_values_for_3_subsets(Data_IA2,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
    print(Check_interim_decision2)
    if (Check_interim_decision2[1] == "continue with 3 subsets") {
      Data_IA3 <- create_data_frame_for_continue_with_3_subsets(q1, q2, q3,
                                                                N_2,
                                                                pi1A,
                                                                pi0A,
                                                                pi1B,
                                                                pi0B, pi1C, pi0C,
                                                                pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA3 <-
        create_data_frame_if_efficacy_subset1.2_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA3 <-
        create_data_frame_if_efficacy_subset1.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA3 <-
        create_data_frame_if_efficacy_subset2.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
    else if (Check_interim_decision2[1] == "efficacy in the subset 1") {
      Data_IA3 <-
        create_data_frame_if_efficacy_subset1_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA3 <-
        create_data_frame_if_efficacy_subset2_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA3 <-
        create_data_frame_if_efficacy_subset3_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      inits_IA3 <-
        create_list_for_inits_values_for_3_subsets(Data_IA3,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
    print(Check_interim_decision3)
    
    if (Check_interim_decision3[1] == "continue with 3 subsets") {
      Data_IA4 <- create_data_frame_for_continue_with_3_subsets(q1, q2, q3,
                                                                N_2,
                                                                pi1A,
                                                                pi0A,
                                                                pi1B,
                                                                pi0B, pi1C, pi0C,
                                                                pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA4 <-
        create_data_frame_if_efficacy_subset1.2_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      
      
    }
    else if (Check_interim_decision3[1] == "efficacy in the subset 1 & 3") {
      Data_IA4 <-
        create_data_frame_if_efficacy_subset1.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA4 <-
        create_data_frame_if_efficacy_subset2.3_3_subsets(q1, q2, q3,
                                                          N_2,
                                                          pi1A,
                                                          pi0A,
                                                          pi1B,
                                                          pi0B, pi1C, pi0C,
                                                          pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      
      
    }
    else if (Check_interim_decision3[1] == "efficacy in the subset 1") {
      Data_IA4 <-
        create_data_frame_if_efficacy_subset1_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA4 <-
        create_data_frame_if_efficacy_subset2_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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
      Data_IA4 <-
        create_data_frame_if_efficacy_subset3_3_subsets(q1, q2, q3,
                                                        N_2,
                                                        pi1A,
                                                        pi0A,
                                                        pi1B,
                                                        pi0B, pi1C, pi0C,
                                                        pa1, pa2, pa3)
      
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      inits_IA4 <-
        create_list_for_inits_values_for_3_subsets(Data_IA4,
                                                   pi1A,
                                                   pi0A,
                                                   pi1B,
                                                   pi0B,
                                                   pi1C,
                                                   pi0C,
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


 ##### functions for subset 2 with Gail and Simon interaction Test#####
stopping_rules_first_analyse_variation_2_subset_GS <-
  function(P1_0, P1_1, P5, P6, gamma, epsilon) {
    nDecision <- "continue with 2 subsets"
    efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
    efficacy_0 <- ifelse(P1_0 > gamma, 1, 0)
    interaction_quali <- ifelse(P5 > epsilon, 1, 0)
    interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
    interaction_type <-
      ifelse(
        interaction_quali == 1 & interaction_quanti == 1,
        "both interaction",
        ifelse(
          interaction_quali == 1,
          "quali",
          ifelse(interaction_quanti == 1, "quanti", "NO")
        )
      )
    enrichment_1 <-
      ifelse(efficacy_1 == 1 &
               (interaction_quali == 1 | interaction_quanti == 1),
             1,
             0)
    enrichment_0 <-
      ifelse(efficacy_0 == 1 &
               (interaction_quali == 1 | interaction_quanti == 1),
             1,
             0)
    continue <- 0
    if (enrichment_1 == 1) {
      nDecision <- "efficacy in the subset 1"
    } else if (enrichment_0) {
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
        interaction_quali,
        interaction_quanti,
        enrichment_0,
        enrichment_1,
        interaction_type
      )
    )
  }

stopping_rules_other_analyse_variation_2_subset_GS <-
  function(nDecision,
           P1_0,
           P1_1,
           P5,
           P6,
           gamma,
           epsilon,
           interaction_type) {
    continue <- 0
    efficacy_1 <- 0
    efficacy_0 <- 0
    interaction_quali <- 0
    interaction_quanti <- 0
    enrichment_1 <- 0
    enrichment_0 <- 0
    if (nDecision == "efficacy in the subset 1") {
      nDecision <- "efficacy in the subset 1"
      enrichment_1 <- 1
      efficacy_1 <- 1
      interaction_quali <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quali"
        ),
        1,
        0)
      interaction_quanti <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quanti"
        ),
        1,
        0)
    }
    else if (nDecision == "efficacy in the subset 0") {
      nDecision <- "efficacy in the subset 0"
      enrichment_0 <- 1
      efficacy_0 <- 1
      interaction_quali <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quali"
        ),
        1,
        0)
      interaction_quanti <-
        ifelse((
          interaction_type == "both interaction" |
            interaction_type == "quanti"
        ),
        1,
        0)
    }
    else {
      efficacy_1 <- ifelse(P1_1 > gamma, 1, 0)
      efficacy_0 <- ifelse(P1_0 > gamma, 1, 0)
      interaction_quali <- ifelse(P5 > epsilon, 1, 0)
      interaction_quanti <- ifelse(P6 > epsilon, 1, 0)
      enrichment_1 <-
        ifelse(efficacy_1 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      enrichment_0 <-
        ifelse(efficacy_0 == 1 &
                 (interaction_quali == 1 | interaction_quanti == 1),
               1,
               0)
      interaction_type <-
        ifelse(
          interaction_quali == 1 & interaction_quanti == 1,
          "both interaction",
          ifelse(
            interaction_quali == 1,
            "quali",
            ifelse(interaction_quanti == 1, "quanti", "NO")
          )
        )
      if (enrichment_1 == 1) {
        nDecision <- "efficacy in the subset 1"
      }
      else if (enrichment_0 == 1) {
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
        interaction_quali,
        interaction_quanti,
        enrichment_0,
        enrichment_1,
        interaction_type
      )
    )
  }

create_matrix_2_subset_GS <- function(nrow) {
  mat.2 <-
    matrix(nrow = nrow,
           ncol = 54,
           dimnames = list(
             NULL,
             c(
               "gamma",
               "C1",
               "C2",
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
           ))
  
}# create matrix for simulation


#### function for simulation data generation with 2 subset and GS interaction
create_data_frame_for_continue_with_2_subsets <- function(q1,
                                                          N_patients,
                                                          pi1A,
                                                          pi1B,
                                                          pi0A,
                                                          pi0B,
                                                          pa) {
  subset <- rbinom(N_patients, 1, pa) # subsets
  q0 <- (0.5 - (q1 * pa)) / (1 - pa)  #randomization
  group <- NULL
  group[subset == 1] <- rbinom(sum(subset == 1), 1, q1)
  group[subset == 0] <-
    rbinom(sum(subset == 0), 1, q0) # Group traitement vs Placebo
  # Subset biomarkers
  y <- NULL                 # Response
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
}#  create data if decision rules aren't met

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
    rbinom(sum(subset == 0), 1, q0)  #Group traitement vs Placebo
  # Subset biomarkers
  y <- NULL                #  Response
  y[group == 1 &
      subset == 0] <-
    rbinom(sum(group == 1 & subset == 0), 1, pi1B)
  y[group == 0 &
      subset == 0] <-
    rbinom(sum(group == 0 & subset == 0), 1, pi0B)
  return(data.frame(group, subset, y))
}  #create data if subset 0 is effective

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
    rbinom(sum(subset == 1), 1, q1)  #Group traitement vs Placebo
   #Subset biomarkers
  y <- NULL                #  Response
  y[group == 1 &
      subset == 1] <-
    rbinom(sum(group == 1 & subset == 1), 1, pi1A)
  y[group == 0 &
      subset == 1] <-
    rbinom(sum(group == 0 & subset == 1), 1, pi0A)
  return(data.frame(group, subset, y))
} # create data if subset 1 is effective



#### Bayesien model for 2 subset and GS interaction in Jags
binary_bayesian_model_2_subset_GS <- function() {
  nb11 ~ dbin(p11, n11)
  nb10 ~ dbin(p10, n10)
  nb01 ~ dbin(p01, n01)
  nb00 ~ dbin(p00, n00)
  
  p10 ~ dbeta(1, 1)
  p00 ~ dbeta(1, 1)
  p01 ~ dbeta(1, 1)
  p11 ~ dbeta(1, 1)
  
  
  p1.estim <- (p11 * pa) + (p10 * (1 - pa)) # probability in treatment group
  p0.estim <- (p01 * pa) + (p00 * (1 - pa))  # probability in control group
  RR <- p1.estim / p0.estim
  RR.subset1 <- p11 / p01  #effect of treatment in subset 1
  RR.subset0 <- p10 / p00  #effect of treatment in subset 0
  P1_1 <- (RR.subset1 < lambda)
  P1_0 <- (RR.subset0 < lambda)
  biais_0 <- RR.subset0 - (pi1B / pi0B)
  biais_1 <- RR.subset1 - (pi1A / pi0A)
  logRR1 <- log(RR.subset1)
  logRR0 <- log(RR.subset0)
  sum1 <- (RR.subset1 - RR)
  sum0 <- (RR.subset0 - RR)
  sum_final <- sum(sum1, sum0)
  P6 <- (sum_final > C2)
  sd1 <- sqrt(mean((y1 - mean(y1)) ^ 2))
  sd0 <- sqrt(mean((y0 - mean(y0)) ^ 2))
  betas <- c(logRR1, logRR0)
  sigmas <- c(sd1, sd0)
  over <- betas / sigmas
  oversq <- over ^ 2
  test_quali <-
    min(sum(oversq * (over > 0)), sum(oversq * (over < 0)))
  P5 <- (test_quali > C1)
  
} # bayesian model set 1 for jags function



### function for initialization values for jags wi th 2 subset and GS interaction
create_list_for_inits_values_2_subset_GS <-
  function(Data,
           pi1A,
           pi1B,
           pi0A,
           pi0B,
           pa,
           lambda,
           C1,
           C2) {
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
      pi0B = pi0B,
      pa = pa,
      lambda = lambda,
      C1 = C1,
      C2 = C2,
      y1 = filter(Data, subset == 1)[, 3],
      y0 = filter(Data, subset == 0)[, 3]
    )
  } #list for iniatialization values for jags



### function for result jag model with 2 subsets and GS interaction
results_bayesian_model_2_subset_GS <-
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
    P5 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P5
    P6 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$P6
    biais_0 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_0
    biais_0_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[8, 2]
    biais_1 <- bayesian_jags_fonction_IA1$BUGSoutput$mean$biais_1
    biais_1_sd <-
      bayesian_jags_fonction_IA1$BUGSoutput$summary[9, 2]
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
        P5,
        P6,
        biais_0,
        biais_0_sd,
        biais_1,
        biais_1_sd,
        N_subset1_stage1,
        N_subset0_stage1
      )
    results
  }# results from bayesian model

##### Simulation with 2 subsets and GS interaction#####

simulate_binary_bayesian_trial_2_subset_GS <-
  function(q1 = 0.5,
           N_1 = 200,
           N_2 = 200,
           N_3 = 200,
           N_4 = 200,
           pi1A = 0.3,
           pi1B = 0.4,
           pi0A = 0.3,
           pi0B = 0.4,
           pa = 0.5,
           bayesian.model = binary_bayesian_model_2_subset_GS,
           stopping_rules.IA = stopping_rules_first_analyse_variation_2_subset_GS,
           stopping_rules.TA = stopping_rules_other_analyse_variation_2_subset_GS,
           lambda = 0.8,
           gamma = 0.95,
           epsilon = 0.05,
           C1,
           C2) {
    mat.sub1 <- create_matrix_2_subset_GS(Nsimu)
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
      create_list_for_inits_values_2_subset_GS(Data_IA1, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1,C2)
    
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
        progress.bar = "none",
        jags.seed = "18121995"
      )
    
    results_IA1 <-
      results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA1, Data_IA1)
    
    First_result <- list(Data_IA1, results_IA1)
    
    Check_interim_decision1 <-
      stopping_rules.IA(First_result[[2]][10],
                        First_result[[2]][11],
                        First_result[[2]][12],
                        First_result[[2]][13],
                        gamma,
                        epsilon)
    print(Check_interim_decision1)
    
    if (Check_interim_decision1[1] == "continue with 2 subsets") {
      Data_IA2 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                N_2,
                                                                pi1A,
                                                                pi1B,
                                                                pi0A,
                                                                pi0B,
                                                                pa)
      Data_IA2 <- rbind(First_result[[1]], Data_IA2)
      
      inits_iA2 <-
        create_list_for_inits_values_2_subset_GS(Data_IA2, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
      bayesian_jags_fonction_IA2  <-
        jags(
          data = inits_iA2,
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA2 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA2, Data_IA2)
      
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
        create_list_for_inits_values_2_subset_GS(Data_IA2, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA2 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA2, Data_IA2)
      
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
        create_list_for_inits_values_2_subset_GS(Data_IA2, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA2 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA2, Data_IA2)
      
      Second_result <- list(Data_IA2, results_IA2)
    }
    
    Check_interim_decision2 <-
      stopping_rules.TA(
        Check_interim_decision1[1],
        Second_result[[2]][10],
        Second_result[[2]][11],
        Second_result[[2]][12],
        Second_result[[2]][13],
        gamma,
        epsilon,
        Check_interim_decision1[9]
      )
    print(Check_interim_decision2)
    if (Check_interim_decision2[1] == "continue with 2 subsets") {
      Data_IA3 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                N_3,
                                                                pi1A,
                                                                pi1B,
                                                                pi0A,
                                                                pi0B,
                                                                pa)
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      
      inits_IA3 <-
        create_list_for_inits_values_2_subset_GS(Data_IA3, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA3 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA3, Data_IA3)
      
      Third_result <- list(Data_IA3, results_IA3)
      
    } else if ((Check_interim_decision2[1] == "efficacy in the subset 0")) {
      Data_IA3 <- create_data_frame_if_efficacy_subset0(q1,
                                                        N_3,
                                                        pi1A,
                                                        pi1B,
                                                        pi0A,
                                                        pi0B,
                                                        pa)
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      
      
      
      inits_IA3 <-
        create_list_for_inits_values_2_subset_GS(Data_IA3, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA3 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA3, Data_IA3)
      
      Third_result <- list(Data_IA3, results_IA3)
      
    } else if ((Check_interim_decision2[1] == "efficacy in the subset 1")) {
      Data_IA3 <- create_data_frame_if_efficacy_subset1(q1,
                                                        N_3,
                                                        pi1A,
                                                        pi1B,
                                                        pi0A,
                                                        pi0B,
                                                        pa)
      Data_IA3 <- rbind(Second_result[[1]], Data_IA3)
      
      inits_IA3 <-
        create_list_for_inits_values_2_subset_GS(Data_IA3, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA3 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA3, Data_IA3)
      
      Third_result <- list(Data_IA3, results_IA3)
    }
    
    Check_interim_decision3 <-
      stopping_rules.TA(
        Check_interim_decision2[1],
        Third_result[[2]][10],
        Third_result[[2]][11],
        Third_result[[2]][12],
        Third_result[[2]][13],
        gamma,
        epsilon,
        Check_interim_decision2[9]
      )
    print(Check_interim_decision3)
    if (Check_interim_decision3[1] == "continue with 2 subsets") {
      Data_IA4 <- create_data_frame_for_continue_with_2_subsets(q1,
                                                                N_4,
                                                                pi1A,
                                                                pi1B,
                                                                pi0A,
                                                                pi0B,
                                                                pa)
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      
      inits_IA4 <-
        create_list_for_inits_values_2_subset_GS(Data_IA4, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA4 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA4, Data_IA4)
      
      Final_result <- list(Data_IA4, results_IA4)
      
    } else if ((Check_interim_decision3[1] == "efficacy in the subset 0")) {
      Data_IA4 <- create_data_frame_if_efficacy_subset0(q1,
                                                        N_4,
                                                        pi1A,
                                                        pi1B,
                                                        pi0A,
                                                        pi0B,
                                                        pa)
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      
      inits_IA4 <-
        create_list_for_inits_values_2_subset_GS(Data_IA4, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA4 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA4, Data_IA4)
      
      Final_result <- list(Data_IA4, results_IA4)
      
    } else if ((Check_interim_decision3[1] == "efficacy in the subset 1")) {
      Data_IA4 <- create_data_frame_if_efficacy_subset1(q1,
                                                        N_4,
                                                        pi1A,
                                                        pi1B,
                                                        pi0A,
                                                        pi0B,
                                                        pa)
      Data_IA4 <- rbind(Third_result[[1]], Data_IA4)
      
      inits_IA4 <-
        create_list_for_inits_values_2_subset_GS(Data_IA4, pi1A, pi1B, pi0A, pi0B, pa, lambda, C1, C2)
      
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
          progress.bar = "none",
          jags.seed = "18121995"
        )
      
      results_IA4 <-
        results_bayesian_model_2_subset_GS(bayesian_jags_fonction_IA4, Data_IA4)
      
      Final_result <- list(Data_IA4, results_IA4)
      
      
    }
    
    Check_interim_decision4 <-
      stopping_rules.TA(
        Check_interim_decision3[1],
        Final_result[[2]][10],
        Final_result[[2]][11],
        Final_result[[2]][12],
        Final_result[[2]][13],
        gamma,
        epsilon,
        Check_interim_decision3[9]
      )
    print(Check_interim_decision4)
    
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
    #}
    return(apply(mat.sub1, 2, function(x) {
      mean(x, na.rm = T)
    }))
  }



#### simulation scenarios #####

#### values for parallelize the scenarios
Nsimu <- 10000
seeds <- sample(1:200000, Nsimu)
gamma_var <- c(0.85,0.90, 0.95)
pa_choice <- c(0.2, 0.4, 0.5, 0.6, 0.8)
q1_choice <- c(0.1, 0.3, 0.5, 0.6, 0.9)
cl <- makeCluster(10,outfile="")
registerDoSNOW(cl)

#### Scenario Simulation with 3 subsets ####
##### scenario1####
#pattern 1
  sc1_3subset_variation_1 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
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
        1 / 3,
        1 / 3,
        1 / 3,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc1_3subset_variation_1[[m]] <- data.frame(res)
  }
  saveRDS(sc1_3subset_variation_1, file = "sc1_3subset_variation_estimate_1.rds")
  sc1_3subset_variation_1 <- bind_rows(sc1_3subset_variation_1)
  
  sc1_3subset_variation_1 <-
    sc1_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc1_3subset_variation_1) <- c(
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
  #pattern 2
  sc1_3subset_variation_2 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
  
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
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
        1 / 6,
        1 / 3,
        2 / 4,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
  
      )
    }
    sc1_3subset_variation_2[[m]] <- data.frame(res)
  }
  saveRDS(sc1_3subset_variation_2, file = "sc1_3subset_variation_estimate_2.rds")
  sc1_3subset_variation_2 <- bind_rows(sc1_3subset_variation_2)
  
  sc1_3subset_variation_2 <-
    sc1_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc1_3subset_variation_2) <- c(
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
 #pattern 3
 sc1_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
 
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
   simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
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
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
 
     )
   }
   sc1_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc1_3subset_variation_3, file = "sc1_3subset_variation_estimate_3.rds")
 sc1_3subset_variation_3 <- bind_rows(sc1_3subset_variation_3)
 
 sc1_3subset_variation_3 <-
   sc1_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc1_3subset_variation_3) <- c(
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
 #pattern 4
  sc1_3subset_variation_4 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
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
        1 / 3,
        2 / 4,
        1 / 6,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc1_3subset_variation_4[[m]] <- data.frame(res)
  }
  saveRDS(sc1_3subset_variation_4, file = "sc1_3subset_variation_estimate_4.rds")
  sc1_3subset_variation_4 <- bind_rows(sc1_3subset_variation_4)
  
  sc1_3subset_variation_4 <-
    sc1_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc1_3subset_variation_4) <- c(
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
 #pattern 5
  sc1_3subset_variation_5 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
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
        1 / 3,
        1 / 18,
        11 / 18,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc1_3subset_variation_5[[m]] <- data.frame(res)
  }
  saveRDS(sc1_3subset_variation_5, file = "sc1_3subset_variation_estimate_5.rds")
  sc1_3subset_variation_5 <- bind_rows(sc1_3subset_variation_5)
  
  sc1_3subset_variation_5 <-
    sc1_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc1_3subset_variation_5) <- c(
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
 
##### scenario 2#####
  #pattern 1
  sc2_3subset_variation_1 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
        200,
        200,
        200,
        200,
        0.32,
        0.4,
        0.32,
        0.4,
        0.32,
        0.4,
        1 / 3,
        1 / 3,
        1 / 3,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc2_3subset_variation_1[[m]] <- data.frame(res)
  }
  saveRDS(sc2_3subset_variation_1, file = "sc2_3subset_variation_estimate_1.rds")
  sc2_3subset_variation_1 <- bind_rows(sc2_3subset_variation_1)
  
  sc2_3subset_variation_1 <-
    sc2_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc2_3subset_variation_1) <- c(
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
  
  #pattern 2
  sc2_3subset_variation_2 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
        200,
        200,
        200,
        200,
        0.32,
        0.4,
        0.32,
        0.4,
        0.32,
        0.4,
        1 / 6,
        1 / 3,
        2 / 4,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc2_3subset_variation_2[[m]] <- data.frame(res)
  }
  saveRDS(sc2_3subset_variation_2, file = "sc2_3subset_variation_estimate_2.rds")
  sc2_3subset_variation_2 <- bind_rows(sc2_3subset_variation_2)
  
  sc2_3subset_variation_2 <-
    sc2_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc2_3subset_variation_2) <- c(
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
 
  #pattern 3
 sc2_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.32,
       0.4,
       0.32,
       0.4,
       0.32,
       0.4,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc2_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc2_3subset_variation_3, file = "sc2_3subset_variation_estimate_3.rds")
 sc2_3subset_variation_3 <- bind_rows(sc2_3subset_variation_3)
 
 sc2_3subset_variation_3 <-
   sc2_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc2_3subset_variation_3) <- c(
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
 
 #pattern 4
  sc2_3subset_variation_4 <- list()
  for (m in seq_along(gamma_var)) {
    gamma_vars <- gamma_var[m]
    
    res <- foreach(
      i = 1:Nsimu,
      .inorder = F,
      .combine = "rbind",
      .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
    ) %dopar% {
      set.seed(seeds[i])
      simulate_binary_bayesian_trial_3_subsets(
        0.5,
        0.5,
        0.5,
        200,
        200,
        200,
        200,
        0.32,
        0.4,
        0.32,
        0.4,
        0.32,
        0.4,
        1 / 3,
        2 / 4,
        1 / 6,
        binary_bayesian_model_3_subset,
        stopping_rules_first_analyse_variation_3_subset,
        stopping_rules_other_analyse_variation_3_subset,
        0.9,
        gamma_vars,
        0.05,
        1.970,
        0.715
        
      )
    }
    sc2_3subset_variation_4[[m]] <- data.frame(res)
  }
  saveRDS(sc2_3subset_variation_4, file = "sc2_3subset_variation_estimate_4.rds")
  sc2_3subset_variation_4 <- bind_rows(sc2_3subset_variation_4)
  
  sc2_3subset_variation_4 <-
    sc2_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
  colnames(sc2_3subset_variation_4) <- c(
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
 
  #pattern 5
 sc2_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.32,
       0.4,
       0.32,
       0.4,
       0.32,
       0.4,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc2_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc2_3subset_variation_5, file = "sc2_3subset_variation_estimate_5.rds")
 sc2_3subset_variation_5 <- bind_rows(sc2_3subset_variation_5)
 
 sc2_3subset_variation_5 <-
   sc2_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc2_3subset_variation_5) <- c(
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
 
 ####scenario 3####
 #pattern 1
 sc3_3subset_variation_1 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.2,
       0.5,
       1 / 3,
       1 / 3,
       1 / 3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc3_3subset_variation_1[[m]] <- data.frame(res)
 }
 saveRDS(sc3_3subset_variation_1, file = "sc3_3subset_variation_estimate_1.rds")
 sc3_3subset_variation_1 <- bind_rows(sc3_3subset_variation_1)
 
 sc3_3subset_variation_1 <-
   sc3_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc3_3subset_variation_1) <- c(
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
 
 #pattern 2
 sc3_3subset_variation_2 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.2,
       0.5,
       1 / 6,
       1 / 3,
       2 / 4,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc3_3subset_variation_2[[m]] <- data.frame(res)
 }
 saveRDS(sc3_3subset_variation_2, file = "sc3_3subset_variation_estimate_2.rds")
 sc3_3subset_variation_2 <- bind_rows(sc3_3subset_variation_2)
 
 sc3_3subset_variation_2 <-
   sc3_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc3_3subset_variation_2) <- c(
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
 #pattern 3
 sc3_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.2,
       0.5,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc3_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc3_3subset_variation_3, file = "sc3_3subset_variation_estimate_3.rds")
 sc3_3subset_variation_3 <- bind_rows(sc3_3subset_variation_3)
 
 sc3_3subset_variation_3 <-
   sc3_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc3_3subset_variation_3) <- c(
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
 #pattern 4
 sc3_3subset_variation_4 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.2,
       0.5,
       1 / 3,
       2 / 4,
       1 / 6,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc3_3subset_variation_4[[m]] <- data.frame(res)
 }
 saveRDS(sc3_3subset_variation_4, file = "sc3_3subset_variation_estimate_4.rds")
 sc3_3subset_variation_4 <- bind_rows(sc3_3subset_variation_4)
 
 sc3_3subset_variation_4 <-
   sc3_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc3_3subset_variation_4) <- c(
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
 #pattern 5
 sc3_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.4,
       0.4,
       0.2,
       0.5,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc3_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc3_3subset_variation_5, file = "sc3_3subset_variation_estimate_5.rds")
 sc3_3subset_variation_5 <- bind_rows(sc3_3subset_variation_5)
 
 sc3_3subset_variation_5 <-
   sc3_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc3_3subset_variation_5) <- c(
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
#### scenario4#####
 #pattern 1
 sc4_3subset_variation_1 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.24,
       0.4,
       1 / 3,
       1 / 3,
       1 / 3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc4_3subset_variation_1[[m]] <- data.frame(res)
 }
 saveRDS(sc4_3subset_variation_1, file = "sc4_3subset_variation_estimate_1.rds")
 sc4_3subset_variation_1 <- bind_rows(sc4_3subset_variation_1)
 
 sc4_3subset_variation_1 <-
   sc4_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc4_3subset_variation_1) <- c(
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
 #pattern 2
 sc4_3subset_variation_2 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.24,
       0.4,
       1 / 6,
       1 / 3,
       2 / 4,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc4_3subset_variation_2[[m]] <- data.frame(res)
 }
 saveRDS(sc4_3subset_variation_2, file = "sc4_3subset_variation_estimate_2.rds")
 sc4_3subset_variation_2 <- bind_rows(sc4_3subset_variation_2)
 
 sc4_3subset_variation_2 <-
   sc4_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc4_3subset_variation_2) <- c(
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
 #pattern 3
 sc4_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.24,
       0.4,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc4_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc4_3subset_variation_3, file = "sc4_3subset_variation_estimate_3.rds")
 sc4_3subset_variation_3 <- bind_rows(sc4_3subset_variation_3)
 
 sc4_3subset_variation_3 <-
   sc4_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc4_3subset_variation_3) <- c(
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
 #pattern 4
 sc4_3subset_variation_4 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.24,
       0.4,
       1 / 3,
       2 / 4,
       1 / 6,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc4_3subset_variation_4[[m]] <- data.frame(res)
 }
 saveRDS(sc4_3subset_variation_4, file = "sc4_3subset_variation_estimate_4.rds")
 sc4_3subset_variation_4 <- bind_rows(sc4_3subset_variation_4)
 
 sc4_3subset_variation_4 <-
   sc4_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc4_3subset_variation_4) <- c(
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
 #pattern 5
 sc4_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.24,
       0.4,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc4_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc4_3subset_variation_5, file = "sc4_3subset_variation_estimate_5.rds")
 sc4_3subset_variation_5 <- bind_rows(sc4_3subset_variation_5)
 
 sc4_3subset_variation_5 <-
   sc4_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc4_3subset_variation_5) <- c(
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
 ####scenario5#####
 #pattern 1
 sc5_3subset_variation_1 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.2,
       0.5,
       1 / 3,
       1 / 3,
       1 / 3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc5_3subset_variation_1[[m]] <- data.frame(res)
 }
 saveRDS(sc5_3subset_variation_1, file = "sc5_3subset_variation_estimate_1.rds")
 sc5_3subset_variation_1 <- bind_rows(sc5_3subset_variation_1)
 
 sc5_3subset_variation_1 <-
   sc5_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc5_3subset_variation_1) <- c(
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
 #pattern 2
 sc5_3subset_variation_2 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.2,
       0.5,
       1 / 6,
       1 / 3,
       2 / 4,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc5_3subset_variation_2[[m]] <- data.frame(res)
 }
 saveRDS(sc5_3subset_variation_2, file = "sc5_3subset_variation_estimate_2.rds")
 sc5_3subset_variation_2 <- bind_rows(sc5_3subset_variation_2)
 
 sc5_3subset_variation_2 <-
   sc5_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc5_3subset_variation_2) <- c(
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
 #pattern 3
 sc5_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.2,
       0.5,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc5_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc5_3subset_variation_3, file = "sc5_3subset_variation_estimate_3.rds")
 sc5_3subset_variation_3 <- bind_rows(sc5_3subset_variation_3)
 
 sc5_3subset_variation_3 <-
   sc5_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc5_3subset_variation_3) <- c(
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
 #pattern 4
 sc5_3subset_variation_4 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.2,
       0.5,
       1 / 3,
       2 / 4,
       1 / 6,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc5_3subset_variation_4[[m]] <- data.frame(res)
 }
 saveRDS(sc5_3subset_variation_4, file = "sc5_3subset_variation_estimate_4.rds")
 sc5_3subset_variation_4 <- bind_rows(sc5_3subset_variation_4)
 
 sc5_3subset_variation_4 <-
   sc5_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc5_3subset_variation_4) <- c(
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
 #pattern 5
 sc5_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.2,
       0.5,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc5_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc5_3subset_variation_5, file = "sc5_3subset_variation_estimate_5.rds")
 sc5_3subset_variation_5 <- bind_rows(sc5_3subset_variation_5)
 
 sc5_3subset_variation_5 <-
   sc5_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc5_3subset_variation_5) <- c(
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

####scenario6#####
 #pattern 1
 sc6_3subset_variation_1 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.5,
       0.4,
       1 / 3,
       1 / 3,
       1 / 3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc6_3subset_variation_1[[m]] <- data.frame(res)
 }
 saveRDS(sc6_3subset_variation_1, file = "sc6_3subset_variation_estimate_1.rds")
 sc6_3subset_variation_1 <- bind_rows(sc6_3subset_variation_1)
 
 sc6_3subset_variation_1 <-
   sc6_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc6_3subset_variation_1) <- c(
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
 #pattern 2
 sc6_3subset_variation_2 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.5,
       0.4,
       1 / 6,
       1 / 3,
       2 / 4,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc6_3subset_variation_2[[m]] <- data.frame(res)
 }
 saveRDS(sc6_3subset_variation_2, file = "sc6_3subset_variation_estimate_2.rds")
 sc6_3subset_variation_2 <- bind_rows(sc6_3subset_variation_2)
 
 sc6_3subset_variation_2 <-
   sc6_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc6_3subset_variation_2) <- c(
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
 #pattern 3
 sc6_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.5,
       0.4,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc6_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc6_3subset_variation_3, file = "sc6_3subset_variation_estimate_3.rds")
 sc6_3subset_variation_3 <- bind_rows(sc6_3subset_variation_3)
 
 sc6_3subset_variation_3 <-
   sc6_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc6_3subset_variation_3) <- c(
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
 #pattern 4
 sc6_3subset_variation_4 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.5,
       0.4,
       1 / 3,
       2 / 4,
       1 / 6,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc6_3subset_variation_4[[m]] <- data.frame(res)
 }
 saveRDS(sc6_3subset_variation_4, file = "sc6_3subset_variation_estimate_4.rds")
 sc6_3subset_variation_4 <- bind_rows(sc6_3subset_variation_4)
 
 sc6_3subset_variation_4 <-
   sc6_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc6_3subset_variation_4) <- c(
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
 #pattern 5
 sc6_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.24,
       0.4,
       0.5,
       0.4,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc6_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc6_3subset_variation_5, file = "sc6_3subset_variation_estimate_5.rds")
 sc6_3subset_variation_5 <- bind_rows(sc6_3subset_variation_5)
 
 sc6_3subset_variation_5 <-
   sc6_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc6_3subset_variation_5) <- c(
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
 
#### scenario7#####
 #pattern 1
 sc7_3subset_variation_1 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.2,
       0.5,
       0.5,
       0.4,
       1 / 3,
       1 / 3,
       1 / 3,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc7_3subset_variation_1[[m]] <- data.frame(res)
 }
 saveRDS(sc7_3subset_variation_1, file = "sc7_3subset_variation_estimate_1.rds")
 sc7_3subset_variation_1 <- bind_rows(sc7_3subset_variation_1)
 
 sc7_3subset_variation_1 <-
   sc7_3subset_variation_1 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc7_3subset_variation_1) <- c(
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
 #pattern 2
 sc7_3subset_variation_2 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.2,
       0.5,
       0.5,
       0.4,
       1 / 6,
       1 / 3,
       2 / 4,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc7_3subset_variation_2[[m]] <- data.frame(res)
 }
 saveRDS(sc7_3subset_variation_2, file = "sc7_3subset_variation_estimate_2.rds")
 sc7_3subset_variation_2 <- bind_rows(sc7_3subset_variation_2)
 
 sc7_3subset_variation_2 <-
   sc7_3subset_variation_2 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc7_3subset_variation_2) <- c(
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
 #pattern 3
 sc7_3subset_variation_3 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.2,
       0.5,
       0.5,
       0.4,
       11 / 18,
       1 / 3,
       1 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc7_3subset_variation_3[[m]] <- data.frame(res)
 }
 saveRDS(sc7_3subset_variation_3, file = "sc7_3subset_variation_estimate_3.rds")
 sc7_3subset_variation_3 <- bind_rows(sc7_3subset_variation_3)
 
 sc7_3subset_variation_3 <-
   sc7_3subset_variation_3 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc7_3subset_variation_3) <- c(
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
 #pattern 4
 sc7_3subset_variation_4 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.2,
       0.5,
       0.5,
       0.4,
       1 / 3,
       2 / 4,
       1 / 6,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc7_3subset_variation_4[[m]] <- data.frame(res)
 }
 saveRDS(sc7_3subset_variation_4, file = "sc7_3subset_variation_estimate_4.rds")
 sc7_3subset_variation_4 <- bind_rows(sc7_3subset_variation_4)
 
 sc7_3subset_variation_4 <-
   sc7_3subset_variation_4 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc7_3subset_variation_4) <- c(
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
 #pattern 5
 sc7_3subset_variation_5 <- list()
 for (m in seq_along(gamma_var)) {
   gamma_vars <- gamma_var[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_3_subsets(
       0.5,
       0.5,
       0.5,
       200,
       200,
       200,
       200,
       0.4,
       0.4,
       0.2,
       0.5,
       0.5,
       0.4,
       1 / 3,
       1 / 18,
       11 / 18,
       binary_bayesian_model_3_subset,
       stopping_rules_first_analyse_variation_3_subset,
       stopping_rules_other_analyse_variation_3_subset,
       0.9,
       gamma_vars,
       0.05,
       1.970,
       0.715
       
     )
   }
   sc7_3subset_variation_5[[m]] <- data.frame(res)
 }
 saveRDS(sc7_3subset_variation_5, file = "sc7_3subset_variation_estimate_5.rds")
 sc7_3subset_variation_5 <- bind_rows(sc7_3subset_variation_5)
 
 sc7_3subset_variation_5 <-
   sc7_3subset_variation_5 %>% group_by(gamma) %>% summarise_all(mean)
 colnames(sc7_3subset_variation_5) <- c(
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
#####scenario 1-4 subset2 and GS interaction#####
 
 ### pa correspond to prevalence subset (sensitivity simulation)
 sc1_2subset_variation_pa <- list()
 for (m in seq_along(pa_choice)) {
   pa_choices <- pa_choice[m]
 
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
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
       0.4,
       pa_choices,
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.655,
       0.220
 
 
     )
   }
   sc1_2subset_variation_pa[[m]] <- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_pa, file = "sc1_2subset_variation_estimate_pa.rds")
 sc1_2subset_variation_pa <- bind_rows(sc1_2subset_variation_pa)
 
 sc1_2subset_variation_pa <-
   sc1_2subset_variation_pa %>% group_by(pa) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_pa) <- c(
   "gamma",
   "C1",
   "C2",
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

sc2_2subset_variation_pa <- list()
for (m in seq_along(pa_choice)) {
  pa_choices <- pa_choice[m]

  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2_subset_GS(
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
      binary_bayesian_model_2_subset_GS,
      stopping_rules_first_analyse_variation_2_subset_GS,
      stopping_rules_other_analyse_variation_2_subset_GS,
      0.9,
      0.9,
      0.05,
      0.655,
      0.220


    )
  }
  sc2_2subset_variation_pa[[m]] <- data.frame(res)
}
saveRDS(sc2_2subset_variation_pa, file = "sc2_2subset_variation_estimate_pa.rds")
sc2_2subset_variation_pa <- bind_rows(sc2_2subset_variation_pa)

sc2_2subset_variation_pa <-
  sc2_2subset_variation_pa %>% group_by(pa) %>% summarise_all(mean)
colnames(sc2_2subset_variation_pa) <- c(
  "gamma",
  "C1",
  "C2",
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

 sc3_2subset_variation_pa <- list()
 for (m in seq_along(pa_choice)) {
   pa_choices <- pa_choice[m]
 
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
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
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.655,
       0.220
 
 
     )
   }
   sc3_2subset_variation_pa[[m]] <- data.frame(res)
 }
 saveRDS(sc3_2subset_variation_pa, file = "sc3_2subset_variation_estimate_pa.rds")
 sc3_2subset_variation_pa <- bind_rows(sc3_2subset_variation_pa)
 
 sc3_2subset_variation_pa <-
   sc3_2subset_variation_pa %>% group_by(pa) %>% summarise_all(mean)
 colnames(sc3_2subset_variation_pa) <- c(
   "gamma",
   "C1",
   "C2",
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
 
 sc4_2subset_variation_pa <- list()
 for (m in seq_along(pa_choice)) {
   pa_choices <- pa_choice[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
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
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.655,
       0.220
       
       
     )
   }
   sc4_2subset_variation_pa[[m]] <- data.frame(res)
 }
 saveRDS(sc4_2subset_variation_pa, file = "sc4_2subset_variation_estimate_pa.rds")
 sc4_2subset_variation_pa <- bind_rows(sc4_2subset_variation_pa)
 
 sc4_2subset_variation_pa <-
   sc4_2subset_variation_pa %>% group_by(pa) %>% summarise_all(mean)
 colnames(sc4_2subset_variation_pa) <- c(
   "gamma",
   "C1",
   "C2",
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
 
 
 #### Q1 correspond to balance of randomization simulation ( sensitivity analysis)
 sc1_2subset_variation_q1 <- list()
 for (m in seq_along(q1_choice)) {
   q1_choices <- q1_choice[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
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
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.655,
       0.220
       
       
     )
   }
   sc1_2subset_variation_q1[[m]] <- data.frame(res)
 }
 saveRDS(sc1_2subset_variation_q1, file = "sc1_2subset_variation_estimate_q1.rds")
 sc1_2subset_variation_q1 <- bind_rows(sc1_2subset_variation_q1)
 
 sc1_2subset_variation_q1 <-
   sc1_2subset_variation_q1 %>% group_by(q1) %>% summarise_all(mean)
 colnames(sc1_2subset_variation_q1) <- c(
   "gamma",
   "C1",
   "C2",
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

sc2_2subset_variation_q1 <- list()
for (m in seq_along(q1_choice)) {
  q1_choices <- q1_choice[m]
  
  res <- foreach(
    i = 1:Nsimu,
    .inorder = F,
    .combine = "rbind",
    .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
  ) %dopar% {
    set.seed(seeds[i])
    simulate_binary_bayesian_trial_2_subset_GS(
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
      binary_bayesian_model_2_subset_GS,
      stopping_rules_first_analyse_variation_2_subset_GS,
      stopping_rules_other_analyse_variation_2_subset_GS,
      0.9,
      0.9,
      0.05,
      0.655,
      0.220
      
      
    )
  }
  sc2_2subset_variation_q1[[m]] <- data.frame(res)
}
saveRDS(sc2_2subset_variation_q1, file = "sc2_2subset_variation_estimate_q1.rds")
sc2_2subset_variation_q1 <- bind_rows(sc2_2subset_variation_q1)

sc2_2subset_variation_q1 <-
  sc2_2subset_variation_q1 %>% group_by(q1) %>% summarise_all(mean)
colnames(sc2_2subset_variation_q1) <- c(
  "gamma",
  "C1",
  "C2",
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

 sc3_2subset_variation_q1 <- list()
 for (m in seq_along(q1_choice)) {
   q1_choices <- q1_choice[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
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
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.665,
       0.220
       
       
     )
   }
   sc3_2subset_variation_q1[[m]] <- data.frame(res)
 }
 saveRDS(sc3_2subset_variation_q1, file = "sc3_2subset_variation_estimate_q1.rds")
 sc3_2subset_variation_q1 <- bind_rows(sc3_2subset_variation_q1)
 
 sc3_2subset_variation_q1 <-
   sc3_2subset_variation_q1 %>% group_by(q1) %>% summarise_all(mean)
 colnames(sc3_2subset_variation_q1) <- c(
   "gamma",
   "C1",
   "C2",
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
 
 sc4_2subset_variation_q1 <- list()
 for (m in seq_along(q1_choice)) {
   q1_choices <- q1_choice[m]
   
   res <- foreach(
     i = 1:Nsimu,
     .inorder = F,
     .combine = "rbind",
     .packages = c("dplyr", "R2jags"),.errorhandling = "remove"
   ) %dopar% {
     set.seed(seeds[i])
     simulate_binary_bayesian_trial_2_subset_GS(
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
       binary_bayesian_model_2_subset_GS,
       stopping_rules_first_analyse_variation_2_subset_GS,
       stopping_rules_other_analyse_variation_2_subset_GS,
       0.9,
       0.9,
       0.05,
       0.655,
       0.220
       
       
     )
   }
   sc4_2subset_variation_q1[[m]] <- data.frame(res)
 }
 saveRDS(sc4_2subset_variation_q1, file = "sc4_2subset_variation_estimate_q1.rds")
 sc4_2subset_variation_q1 <- bind_rows(sc4_2subset_variation_q1)
 
 sc4_2subset_variation_q1 <-
   sc4_2subset_variation_q1 %>% group_by(q1) %>% summarise_all(mean)
 colnames(sc4_2subset_variation_q1) <- c(
   "gamma",
   "C1",
   "C2",
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

#### save all scenarios for create graphics and tables #####
 saveRDS(sc1_3subset_variation_1, file = "sc1_3subset_variation_final_1.rds")
 saveRDS(sc1_3subset_variation_2, file = "sc1_3subset_variation_final_2.rds")
 saveRDS(sc1_3subset_variation_3, file = "sc1_3subset_variation_final_3.rds")
 saveRDS(sc1_3subset_variation_4, file = "sc1_3subset_variation_final_4.rds")
 saveRDS(sc1_3subset_variation_5, file = "sc1_3subset_variation_final_5.rds")

 saveRDS(sc2_3subset_variation_1, file = "sc2_3subset_variation_final_1.rds")
 saveRDS(sc2_3subset_variation_2, file = "sc2_3subset_variation_final_2.rds")
 saveRDS(sc2_3subset_variation_3, file = "sc2_3subset_variation_final_3.rds")
 saveRDS(sc2_3subset_variation_4, file = "sc2_3subset_variation_final_4.rds")
 saveRDS(sc2_3subset_variation_5, file = "sc2_3subset_variation_final_5.rds")

 saveRDS(sc3_3subset_variation_1, file = "sc3_3subset_variation_final_1.rds")
 saveRDS(sc3_3subset_variation_2, file = "sc3_3subset_variation_final_2.rds")
 saveRDS(sc3_3subset_variation_3, file = "sc3_3subset_variation_final_3.rds")
 saveRDS(sc3_3subset_variation_4, file = "sc3_3subset_variation_final_4.rds")
 saveRDS(sc3_3subset_variation_5, file = "sc3_3subset_variation_final_5.rds")
 
 saveRDS(sc4_3subset_variation_1, file = "sc4_3subset_variation_final_1.rds")
 saveRDS(sc4_3subset_variation_2, file = "sc4_3subset_variation_final_2.rds")
 saveRDS(sc4_3subset_variation_3, file = "sc4_3subset_variation_final_3.rds")
 saveRDS(sc4_3subset_variation_4, file = "sc4_3subset_variation_final_4.rds")
 saveRDS(sc4_3subset_variation_5, file = "sc4_3subset_variation_final_5.rds")
 
 saveRDS(sc5_3subset_variation_1, file = "sc5_3subset_variation_final_1.rds")
 saveRDS(sc5_3subset_variation_2, file = "sc5_3subset_variation_final_2.rds")
 saveRDS(sc5_3subset_variation_3, file = "sc5_3subset_variation_final_3.rds")
 saveRDS(sc5_3subset_variation_4, file = "sc5_3subset_variation_final_4.rds")
 saveRDS(sc5_3subset_variation_5, file = "sc5_3subset_variation_final_5.rds")
 
 saveRDS(sc6_3subset_variation_1, file = "sc6_3subset_variation_final_1.rds")
 saveRDS(sc6_3subset_variation_2, file = "sc6_3subset_variation_final_2.rds")
 saveRDS(sc6_3subset_variation_3, file = "sc6_3subset_variation_final_3.rds")
 saveRDS(sc6_3subset_variation_4, file = "sc6_3subset_variation_final_4.rds")
 saveRDS(sc6_3subset_variation_5, file = "sc6_3subset_variation_final_5.rds")
 
 saveRDS(sc7_3subset_variation_1, file = "sc7_3subset_variation_final_1.rds")
 saveRDS(sc7_3subset_variation_2, file = "sc7_3subset_variation_final_2.rds")
 saveRDS(sc7_3subset_variation_3, file = "sc7_3subset_variation_final_3.rds")
 saveRDS(sc7_3subset_variation_4, file = "sc7_3subset_variation_final_4.rds")
 saveRDS(sc7_3subset_variation_5, file = "sc7_3subset_variation_final_5.rds")

saveRDS(sc1_2subset_variation_pa, file = "sc1_2subset_variation_pa.rds")
saveRDS(sc2_2subset_variation_pa, file = "sc2_2subset_variation_pa.rds")
saveRDS(sc3_2subset_variation_pa, file = "sc3_2subset_variation_pa.rds")
saveRDS(sc4_2subset_variation_pa, file = "sc4_2subset_variation_pa.rds")

saveRDS(sc1_2subset_variation_q1, file = "sc1_2subset_variation_q1.rds")
saveRDS(sc2_2subset_variation_q1, file = "sc2_2subset_variation_q1.rds")
saveRDS(sc3_2subset_variation_q1, file = "sc3_2subset_variation_q1.rds")
saveRDS(sc4_2subset_variation_q1, file = "sc4_2subset_variation_q1.rds")

