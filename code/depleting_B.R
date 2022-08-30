#!/usr/bin/env python3

# FILE TO READ IN RESULTS FROM NAIVE_B AND ANALYSE STATISTICALLY

# Load data
df <- read.csv("../results/depleting_nutrient_B.csv", stringsAsFactors = T)
counts <- read.csv("../results/depleting_B_solutions_count.csv", stringsAsFactors = T)

##### SOLUTION COUNTS ####

# 1) Number of solutions across population size
cor.test(counts$N, counts$N_solutions)

# 2) Number of solutions across correlation_length
cor.test(counts$correlation_length, counts$N_solutions)

# 3) Number of solutions across correlation_time
cor.test(counts$correlation_time, counts$N_solutions)

# 4) Linear model
model <- lm(counts$N_solutions~counts$N+counts$correlation_length+counts$correlation_time)
model <- lm(counts$N_solutions~counts$correlation_length*counts$correlation_time)

##### SURVIVAL TIME #######
# 1) Survival time by enviro
cor.test(df$survival_mean, df$N)
cor.test(df$survival_mean, df$correlation_length)
cor.test(df$survival_mean, df$correlation_time)

model <- lm(df$survival_mean~df$N*df$correlation_length*df$correlation_time)

model <- lm(df$survival_mean~df$correlation_length*df$correlation_time)



# 2) Survival time by parameters
cor.test(df$survival_mean, df$R)
cor.test(df$survival_mean, df$T)
cor.test(df$survival_mean, df$m)
cor.test(df$survival_mean, df$s)

# 3) Survival time by forms 
model <- lm(df$survival_mean~ df$qualitative)


###### Parameters by enviro #########
cor.test(df[c(3:6, 15:16,19)], use="pairwise")

cor.test(df$m, df$correlation_length)
cor.test(df$s, df$correlation_length)
cor.test(df$T, df$correlation_time)

####### Forms ##### 
k.mean <- kmeans(df[3:6], centers=10, start=1000)