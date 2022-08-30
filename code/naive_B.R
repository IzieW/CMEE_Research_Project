#!/usr/bin/env python3

# FILE TO READ IN RESULTS FROM NAIVE_B AND ANALYSE STATISTICALLY

# Load data
times <- read.csv("../results/naive_B_means.csv", stringsAsFactors = T)
df <- read.csv("../results/naive_nutrient_B.csv", stringsAsFactors = T)
counts <- read.csv("../results/naive_B_solution_count.csv", stringsAsFactors = T)

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
cor(df[c(2:5, 14:15)], use="pairwise")

cor.test(df$m, df$correlation_length)
cor.test(df$s, df$correlation_length)
cor.test(df$T, df$correlation_time)

####### Forms ##### 
k.mean <- kmeans(df[3:6], centers=10, start=1000)


### do test ###

get_cor_parameters_l <- function(dat=df){
  t <- c(
    cor.test(dat$R, dat$correlation_length)$p.value <0.05,
    cor.test(dat$T, dat$correlation_length)$p.value <0.05,
    cor.test(dat$m, dat$correlation_length)$p.value <0.05,
    cor.test(dat$s, dat$correlation_length)$p.value <0.05
  )
  
  return(t)
}

get_cor_parameters_t <- function(dat=df){
  t <- c(
    cor.test(dat$R, dat$correlation_time)$p.value <0.05,
    cor.test(dat$T, dat$correlation_time)$p.value <0.05,
    cor.test(dat$m, dat$correlation_time)$p.value <0.05,
    cor.test(dat$s, dat$correlation_time)$p.value <0.05
  )
  
  return(t)
}



get_cor <- function(dat=df){
  theta <- c("R", "T", "m", "s")
  d <- data.frame(theta, 
                  get_cor_parameters_l(dat),
                  get_cor_parameters_t(dat))
  colnames(d) <- c("parameter", "correlation_length", "correlation_time")
  
  return(d)
}






#### GROUPS ###
"do_test <- function(dat1, dat2){
  r <- var.test(dat1$R, dat2$R)
  t <- var.test(dat1$T, dat2$T)
  s <- var.test(dat1$s, dat2$s)
  m <- var.test(dat1$m, dat2$m)
  mass <- var.test(dat1$mean_mass, dat2$mean_mass)
  volume <- var.test(dat1$mean_volume, dat2$mean_volume)
  survival <- var.test(dat1$survival_mean, dat2$survival_mean)
  return(c(r$p.value<0.05, t$p.value<0.05, s$p.value<0.05, m$p.value<0.05, survival$p.value<0.05, mass$p.value<0.05, volume$p.value<0.05))
}
"





test_group <- function(dat=df){
  # Ftest to see if groups exist along given parameter
  dat = dat %>% filter(qualitative != "NA")
  forms <- unique(dat$qualitative)
  
  rno = 7 # r number
  tab = matrix(1:8*length(forms), length(forms), 8)
  i=1
  for(form in forms){
    d1 = filter(df, qualitative == form)
    d2 = filter(df, qualitative != form)
    tab[i,] <- append(form, do_test(d1, d2))
    i = i+1
  }
  
  tab = as.data.frame(tab)
  
  colnames(tab) <- c("qualitative", "R", "T", "m", "s", "survival_time", "mass", "volume")
  
  return(tab)
}


with(gap.stat, plot(colMeans(DDwGap),pch=15,type='b',
                    ylim=extendrange(colMeans(DDwGap),f=0.2),
                    xlab="Number of Clusters", ylab="Weighted Gap Statistic"))

