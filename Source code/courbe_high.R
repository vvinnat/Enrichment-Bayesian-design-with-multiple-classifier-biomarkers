#script donnees high
library(tidyverse)
library(R2jags)
library(dplyr)
require(rmeta)

### global values ####
#load Data High
data_HIGH # load Data_high or read.csv function from the repository

####probability of death in each HIGH partition####
####Age partition####
#subset B
#effectif in each treatment group
n0A_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE >= 65, na.rm = T)
n0B_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE < 65, na.rm = T)
#response in each treatment group
r0A_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE >= 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r0B_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE < 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)

# subset  A
#effectif in each treatment group
n1A_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE >= 65, na.rm = T)
n1B_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE < 65, na.rm = T)
#response in each treatment group
r1A_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE >= 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r1B_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE < 65 & data_HIGH$high.dc28 == 1,
      na.rm = T)

pi1A_HIGH_age <-
  r1A_age / n1A_age ## proportion of death in group 1 and subset A
pi1B_HIGH_age <-
  r1B_age / n1B_age ## proportion of death in group 1 and subset B
pi0A_HIGH_age <-
  r0A_age / n0A_age ## proportion of death in group 0 and subset A
pi0B_HIGH_age <-
  r0B_age / n0B_age ## proportion of death in group 0 and subset B
pa_HIGH_age <-
  mean(data_HIGH$high.AGE >= 65) ## prevalence of subset A


###Sofa neuro####
#subset B
n0A_neuro <- sum(data_HIGH$high.rando == 2 &
                   data_HIGH$high.SOFAN >= 1, na.rm = T)
n0B_neuro <- sum(data_HIGH$high.rando == 2 &
                   data_HIGH$high.SOFAN < 1, na.rm = T)
#response in each treatment group
r0A_neuro <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.SOFAN >= 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r0B_neuro <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.SOFAN < 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)

# subset A
#effectif in each treatment group
n1A_neuro <- sum(data_HIGH$high.rando == 1 &
                   data_HIGH$high.SOFAN >= 1, na.rm = T)
n1B_neuro <- sum(data_HIGH$high.rando == 1 &
                   data_HIGH$high.SOFAN < 1, na.rm = T)
#response in each treatment group
r1A_neuro <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.SOFAN >= 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r1B_neuro <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.SOFAN < 1 & data_HIGH$high.dc28 == 1,
      na.rm = T)

pi1A_HIGH_sofa_neuro <- r1A_neuro / n1A_neuro
pi1B_HIGH_sofa_neuro <- r1B_neuro / n1B_neuro
pi0A_HIGH_sofa_neuro <- r0A_neuro / n0A_neuro
pi0B_HIGH_sofa_neuro <- r0B_neuro / n0B_neuro
pa_HIGH_sofa_neuro <- mean(data_HIGH$high.SOFAN >= 1)


###Pao2/fio2 ####
# subset B
n0A_respi <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.sofa_resp == 4, na.rm =
        T)
n0B_respi <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.sofa_resp == 3, na.rm =
        T)
#response in each treatment group
r0A_respi <-
  sum(
    data_HIGH$high.rando == 2 &
      data_HIGH$high.sofa_resp == 4 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )
r0B_respi <-
  sum(
    data_HIGH$high.rando == 2 &
      data_HIGH$high.sofa_resp == 3 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )

# subset A
#effectif in each treatment group
n1A_respi <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.sofa_resp == 4, na.rm =
        T)
n1B_respi <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.sofa_resp == 3, na.rm =
        T)
#response in each treatment group
r1A_respi <-
  sum(
    data_HIGH$high.rando == 1 &
      data_HIGH$high.sofa_resp == 4 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )
r1B_respi <-
  sum(
    data_HIGH$high.rando == 1 &
      data_HIGH$high.sofa_resp == 3 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )

pi1A_HIGH_sofa_resp <- r1A_respi / n1A_respi
pi1B_HIGH_sofa_resp <- r1B_respi / n1B_respi
pi0A_HIGH_sofa_resp <- r0A_respi / n0A_respi
pi0B_HIGH_sofa_resp <- r0B_respi / n0B_respi
pa_HIGH_sofa_resp <- mean(data_HIGH$high.sofa_resp == 4)


### Age with 3 susbets
table(cut(data_HIGH$high.AGE, c(18, 58, 68, 90), include.lowest = TRUE))

#subset B
n1A_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE <= 58, na.rm = T)
n0A_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE <= 58, na.rm = T)
#response in each treatment group

r1A_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE <= 58 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r0A_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE <= 58 & data_HIGH$high.dc28 == 1,
      na.rm = T)

# subset A
#effectif in each treatment group
n1B_age <- sum(data_HIGH$high.rando == 1 & data_HIGH$high.AGE > 58 &
                 data_HIGH$high.AGE <= 68,
               na.rm = T)
n0B_age <- sum(data_HIGH$high.rando == 2 & data_HIGH$high.AGE > 58 &
                 data_HIGH$high.AGE <= 68,
               na.rm = T)
#response in each treatment group
r1B_age <-
  sum(
    data_HIGH$high.rando == 1 & data_HIGH$high.AGE > 58 &
      data_HIGH$high.AGE <= 68 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )
r0B_age <-
  sum(
    data_HIGH$high.rando == 2 & data_HIGH$high.AGE > 58 &
      data_HIGH$high.AGE <= 68 & data_HIGH$high.dc28 == 1,
    na.rm = T
  )

n1C_age <- sum(data_HIGH$high.rando == 1 &
                 data_HIGH$high.AGE > 68, na.rm = T)
n0C_age <- sum(data_HIGH$high.rando == 2 &
                 data_HIGH$high.AGE > 68, na.rm = T)
#response in each treatment group
r1C_age <-
  sum(data_HIGH$high.rando == 1 &
        data_HIGH$high.AGE > 68 & data_HIGH$high.dc28 == 1,
      na.rm = T)
r0C_age <-
  sum(data_HIGH$high.rando == 2 &
        data_HIGH$high.AGE > 68 & data_HIGH$high.dc28 == 1,
      na.rm = T)


pi1A_HIGH_age_3subset <- r1A_age / n1A_age
pi0A_HIGH_age_3subset <- r1B_age / n1B_age
pi1B_HIGH_age_3subset <- r1B_age / n1B_age
pi0B_HIGH_age_3subset <- r0B_age / n0B_age
pi1C_HIGH_age_3subset <- r1C_age / n1C_age
pi0C_HIGH_age_3subset <- r0C_age / n0C_age
pa_HIGH_age_sub1 <- mean(data_HIGH$high.AGE <= 58)
pa_HIGH_age_sub2 <-
  mean(data_HIGH$high.AGE > 58 & data_HIGH$high.AGE <= 68)
pa_HIGH_age_sub3 <- mean(data_HIGH$high.AGE > 68)




####HIGH curve Beta law ######

#### PAO2/FIO2 partition ####
ggplot() + xlim(0, 1) + geom_function(
  aes(linetype = "A", color = "A"),
  size = 1,
  fun = dbeta,
  args = list(shape1 = r1A_respi, shape2 = n1A_respi - r1A_respi)
) +
  geom_function(
    aes(linetype = "B", color = "B"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r1B_respi, shape2 = n1B_respi - r1B_respi)
  ) +
  geom_function(
    aes(linetype = "C", color = "C"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0A_respi, shape2 = n0A_respi - r0A_respi)
  ) +
  geom_function(
    aes(linetype = "D", color = "D"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0B_respi, shape2 = n0B_respi - r0B_respi)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Pr 28 day death") +
  ylab("Density") + ggtitle(expression(PaO[2] / FiO[2] ~ " Partition")) + theme_classic() + scale_x_continuous(breaks  =
                                                                                                                 seq(0, 1, 0.1)) +
  scale_y_continuous(breaks  =
                       seq(0, 20, 1)) +
  theme(
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 18)
  ) + labs(linetype = "", color = "") +
  scale_linetype_manual(
    values = c("solid", "solid", "dashed", "dashed"),
    labels = c(
      expression(PaO[2] / FiO[2] ~ "" < 100 ~ "," ~ "HPNO," ~ "n=131"),
      expression(PaO[2] / FiO[2] ~ "" >= 100 ~ "," ~ "HPNO," ~ "n=257"),
      expression(PaO[2] / FiO[2] ~ "" < 100 ~ "," ~ O[2] ~ "," ~ "n=147"),
      expression(PaO[2] / FiO[2] ~ "" >= 100 ~ "," ~ O[2] ~ "," ~ "n=241")
    )
  ) +
  scale_color_manual(
    values = rep(c("#E64B35FF", "#4DBBD5FF"), 2),
    labels = c(
      expression(PaO[2] / FiO[2] ~ "" < 100 ~ "," ~ "HPNO," ~ "n=131"),
      expression(PaO[2] / FiO[2] ~ "" >= 100 ~ "," ~ "HPNO," ~ "n=257"),
      expression(PaO[2] / FiO[2] ~ "" < 100 ~ "," ~ O[2] ~ "," ~ "n=147"),
      expression(PaO[2] / FiO[2] ~ "" >= 100 ~ "," ~ O[2] ~ "," ~ "n=241")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text.align = 0,
    legend.background = element_rect(
      size = 0.5,
      linetype = "solid",
      colour = "black"
    ),
  )


### Age partition with 2 susbet #####
ggplot() + xlim(0, 1) + geom_function(
  aes(linetype = "A", color = "A"),
  size = 1,
  fun = dbeta,
  args = list(shape1 = r1A_age, shape2 = n1A_age - r1A_age)
) +
  geom_function(
    aes(linetype = "B", color = "B"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r1B_age, shape2 = n1B_age - r1B_age)
  ) +
  geom_function(
    aes(linetype = "C", color = "C"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0A_age, shape2 = n0A_age - r0A_age)
  ) +
  geom_function(
    aes(linetype = "D", color = "D"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0B_age, shape2 = n0B_age - r0B_age)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Pr 28 day death") +
  ylab("Density") + ggtitle("Age Partition with two subets") + theme_classic() + scale_x_continuous(breaks  =
                                                                                                      seq(0, 1, 0.1)) +
  scale_y_continuous(breaks  =
                       seq(0, 20, 1)) +
  theme(
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 18)
  ) + labs(linetype = "", color = "") +
  scale_linetype_manual(
    values = c("solid", "solid", "dashed", "dashed"),
    labels = c(
      expression(Age ~ "" >= 65 ~ "," ~ "HPNO," ~ "n=181"),
      expression(Age ~ "" < 65 ~ "," ~ "HPNO," ~ "n=207"),
      expression(Age ~ "" >= 65 ~ "" ~ "," ~ O[2] ~ "," ~ "n=179"),
      expression(Age ~ "" < 65 ~ "," ~ O[2] ~ "," ~ "n=209")
    )
  ) +
  scale_color_manual(
    values = rep(c("#E64B35FF", "#4DBBD5FF"), 2),
    labels = c(
      expression(Age ~ "" >= 65 ~ "," ~ "HPNO," ~ "n=181"),
      expression(Age ~ "" < 65 ~ "," ~ "HPNO," ~ "n=207"),
      expression(Age ~ "" >= 65 ~ "" ~ "," ~ O[2] ~ "," ~ "n=179"),
      expression(Age ~ "" < 65 ~ "," ~ O[2] ~ "," ~ "n=209")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text.align = 0,
    legend.background = element_rect(
      size = 0.5,
      linetype = "solid",
      colour = "black"
    ),
  )


##### Neurological sofa partition #####
ggplot() + xlim(0, 1) + geom_function(
  aes(linetype = "A", color = "A"),
  size = 1,
  fun = dbeta,
  args = list(shape1 = r1A_neuro, shape2 = n1A_neuro - r1A_neuro)
) +
  geom_function(
    aes(linetype = "B", color = "B"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r1B_neuro, shape2 = n1B_neuro - r1B_neuro)
  ) +
  geom_function(
    aes(linetype = "C", color = "C"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0A_neuro, shape2 = n0A_neuro - r0A_neuro)
  ) +
  geom_function(
    aes(linetype = "D", color = "D"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0B_neuro, shape2 = n0B_neuro - r0B_neuro)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Pr 28 day death") +
  ylab("Density") + ggtitle("SOFA Neurological Partition") + theme_classic() + scale_x_continuous(breaks  =
                                                                                                    seq(0, 1, 0.1)) +
  scale_y_continuous(breaks  =
                       seq(0, 20, 1)) +
  theme(
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 18)
  ) + labs(linetype = "", color = "") +
  scale_linetype_manual(
    values = c("solid", "solid", "dashed", "dashed"),
    labels = c(
      expression("SOFA Neurological" ~ "" >= 1 ~ "," ~ "HPNO," ~ "n=51"),
      expression("SOFA Neurological" ~ "" < 1 ~ "," ~ "HPNO," ~ "n=337"),
      expression("SOFA Neurological" ~ "" >= 1 ~ "," ~ O[2] ~ "," ~ "n=51"),
      expression("SOFA Neurological" ~ "" < 1 ~ "," ~ O[2] ~ "," ~ "n=337")
    )
  ) +
  scale_color_manual(
    values = rep(c("#E64B35FF", "#4DBBD5FF"), 2),
    labels = c(
      expression("SOFA Neurological" ~ "" >= 1 ~ "," ~ "HPNO," ~ "n=51"),
      expression("SOFA Neurological" ~ "" < 1 ~ "," ~ "HPNO," ~ "n=337"),
      expression("SOFA Neurological" ~ "" >= 1 ~ "," ~ O[2] ~ "," ~ "n=51"),
      expression("SOFA Neurological" ~ "" < 1 ~ "," ~ O[2] ~ "," ~ "n=337")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text.align = 0,
    legend.background = element_rect(
      size = 0.5,
      linetype = "solid",
      colour = "black"
    ),
  )

#### Age partition with 3 susbets #####
ggplot() + xlim(0, 1) + geom_function(
  aes(linetype = "A", color = "A"),
  size = 1,
  fun = dbeta,
  args = list(shape1 = r1A_age, shape2 = n1A_age - r1A_age)
) +
  geom_function(
    aes(linetype = "B", color = "B"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r1B_age, shape2 = n1B_age - r1B_age)
  ) +
  geom_function(
    aes(linetype = "C", color = "C"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r1C_age, shape2 = n1C_age - r1C_age)
  ) +
  geom_function(
    aes(linetype = "D", color = "D"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0A_age, shape2 = n0A_age - r0A_age)
  ) +
  geom_function(
    aes(linetype = "E", color = "E"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0B_age, shape2 = n0B_age - r0B_age)
  ) +
  geom_function(
    aes(linetype = "F", color = "F"),
    size = 1,
    fun = dbeta,
    args = list(shape1 = r0C_age, shape2 = n0C_age - r0C_age)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Pr 28 day death") +
  ylab("Density") + ggtitle("Age Partition with three subets") + theme_classic() + scale_x_continuous(breaks  =
                                                                                                        seq(0, 1, 0.1)) +
  scale_y_continuous(breaks  =
                       seq(0, 20, 1)) +
  theme(
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 18)
  ) + labs(linetype = "", color = "") +
  scale_linetype_manual(
    values = c("solid", "solid", "solid", "dashed", "dashed", "dashed"),
    labels = c(
      expression(Age ~ "" <= 58 ~ "," ~ "HPNO," ~ "n=131"),
      expression(Age ~ "" <= 58 ~ "&" > 68 ~ "," ~ "HPNO," ~ "n=129"),
      expression(Age ~ "" > 68 ~ "," ~ "HPNO," ~ "n=128"),
      expression(Age ~ "" <= 58 ~ "" ~ "," ~ O[2] ~ "," ~ "n=124"),
      expression(Age ~ "" <= 58 ~ "&" > 68 ~ "," ~ O[2] ~ "," ~ "n=139"),
      expression(Age ~ "" > 68 ~ "" ~ "," ~ O[2] ~ "," ~ "n=125")
    )
  ) +
  scale_color_manual(
    values = rep(c("#E64B35FF", "#4DBBD5FF", "#00A087FF"), 2),
    labels = c(
      expression(Age ~ "" <= 58 ~ "," ~ "HPNO," ~ "n=131"),
      expression(Age ~ "" <= 58 ~ "&" > 68 ~ "," ~ "HPNO," ~ "n=129"),
      expression(Age ~ "" > 68 ~ "," ~ "HPNO," ~ "n=128"),
      expression(Age ~ "" <= 58 ~ "" ~ "," ~ O[2] ~ "," ~ "n=124"),
      expression(Age ~ "" <= 58 ~ "&" > 68 ~ "," ~ O[2] ~ "," ~ "n=139"),
      expression(Age ~ "" > 68 ~ "" ~ "," ~ O[2] ~ "," ~ "n=125")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text.align = 0,
    legend.background = element_rect(
      size = 0.5,
      linetype = "solid",
      colour = "black"
    ),
    legend.margin = margin(3, 10, 3, 3)
  )
