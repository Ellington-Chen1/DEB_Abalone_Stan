library(deSolve)
library(ggplot2)
library(rstan)

### Simulated Data ###

alpha = exp(-6.36)
beta = 0.87
Resp <- function(Mass){
  R = alpha * Mass ^ beta
  return(R)
}

Mass <- seq(30,80,0.5)

Sim_Data <- Resp(Mass)

gResp <- rgamma(length(Mass), Sim_Data, rate = 400)

Stan_Resp_Data <- list(
  N_obs = length(Mass),
  Mass = Mass,
  Resp = Sim_Data
)

### Model Testing ###

setwd("/Users/homefolder")
Resp_Model_V1 <- stan_model("Metabolism_V1.stan")

test <- sampling(Resp_Model_V1, data = Stan_Resp_Data,
                 iter = 2000,
                 warmup= 1000,
                 chains= 4,
                 cores= 4)

sim_results <- extract(test)

mu <- apply(sim_results$R_hat, 2, mean)

mu_summary <- data.frame(t(apply(sim_results$R_hat, 2, quantile, probs= c(0.025,0.5,0.975))))
names(mu_summary) <- c("CIL","median","CIU")

Comparison_DF <- data.frame(Mass, mu, Sim_Data, mu_summary)

ggplot(Comparison_DF,
       aes(x = Mass)) +
  geom_ribbon(aes(ymin= CIL, ymax= CIU, y= median),alpha= 0.5, fill = "red")+
  geom_line(aes(y = mu), color = "red") +
  geom_point(aes(y = Sim_Data), color = "blue")

hist(sim_results$alpha, breaks = 100)
abline(v = mean(sim_results$alpha), col = "red")
hist(sim_results$log_alpha)
hist(sim_results$beta, breaks = 100)
abline(v = mean(sim_results$beta), col = "red")
hist(sim_results$lambda)

### Real Data ###

setwd("/Users/homefolder/Desktop/Okamoto/DEB_Class")
Mass_OCR <- read.csv("Simple_Mass_OCR.csv")

Real_Stan_Resp_Data <- list(
  N_obs = length(Mass_OCR$Weight),
  Mass = Mass_OCR$Weight,
  Resp = Mass_OCR$Oxygen_Consumption_Rate
)

### Real Model Testing ###

setwd("/Users/homefolder")
Resp_Model_V1 <- stan_model("Metabolism_V1.stan")

Model_V1 <- sampling(Resp_Model_V1, data = Real_Stan_Resp_Data,
                 iter = 4000,
                 warmup= 2000,
                 chains= 4,
                 cores= 4)

results <- extract(Model_V1)

mu <- apply(results$R_hat, 2, mean)

mu_summary <- data.frame(t(apply(results$R_hat, 2, quantile, probs= c(0.025,0.5,0.975))))
names(mu_summary) <- c("CIL","median","CIU")

Comparison_DF <- data.frame(Mass_OCR$Weight, mu, Mass_OCR$Oxygen_Consumption_Rate, mu_summary)

ggplot(Comparison_DF,
       aes(x = Mass_OCR.Weight)) +
  geom_ribbon(aes(ymin= CIL, ymax= CIU, y= median),alpha= 0.5, fill = "red")+
  geom_line(aes(y = mu), color = "red") +
  geom_point(aes(y = Mass_OCR.Oxygen_Consumption_Rate), color = "blue")

hist(results$alpha, breaks = 100)
abline(v = median(results$alpha), col = "red")
hist(results$log_alpha)
hist(results$beta, breaks = 100)
abline(v = mean(results$beta), col = "red")
hist(results$lambda)

