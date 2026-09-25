################################################
### Template for Assignment 3                ###
### MEF 2025-2026                            ###
################################################

rm(list=ls())
library(quadprog)
library(gmm)
library(ggplot2)

#################
### Load Data ###
#################

Dataset     <- read.csv("C:/Users/fvane/Downloads/Assignment_3_dataset.csv")
t              <- dim(Dataset)[1]     # Number of observations
Returns        <- Dataset[, 7:31]     # Rate of returns of the 25 test portfolios
N              <- ncol(Returns)       # N = 25

######################################################################
### Question 1: Expected Utility frontier via quadprog             ###
###                                                                ###
### CARA: min  (lambda/2) w'Sigma w - mu'w   sub w'1 = 1           ###
###                                                                ###
### solve.QP minimizes (1/2) w'D w - d'w    sub A'w=b, so we take  ###
###   D = lambda*Sigma,     d = mu,                                 ###
###   A = 1,                 b = 1,        meq = 1                  ###
######################################################################

### Sample mean vector and covariance matrix of returns
muHat  <- colMeans(Returns)
SigHat <- cov(Returns)

### Grid of lambda: log(lambda) equally spaced on [0.01,4], 250 points
NLam      <- 250
LogLambda <- seq(from=0.01, to=4, length.out=NLam)
Lambda    <- exp(LogLambda)

### Constraint: w'1 = 1  (equality)
A <- matrix(1, nrow = N, ncol = 1)
b <- 1

### Initialize vectors to store portfolio mean and standard deviation
mu_EU <- rep(0, NLam)
sd_EU <- rep(0, NLam)

### For each lambda, solve the optimization problem and store the results
for(i in 1:NLam){
  sol <- solve.QP(Dmat = Lambda[i] * SigHat,
                  dvec = muHat,
                  Amat = A,
                  bvec = b,
                  meq  = 1)
  
  w        <- sol$solution
  mu_EU[i] <- sum(w*muHat)
  sd_EU[i] <- sqrt(as.numeric(t(w) %*% SigHat %*% w))
}

### Plot the EU frontier
plot(sd_EU, mu_EU, type="l", lwd=2,
     xlab="Portfolio Standard Deviation",
     ylab="Portfolio Mean Return",
     main="Expected-Utility Frontier")

dev.copy(png, filename="EU_frontier.png", width=800, height=600)
dev.off()

################################################################
### Question 2: Markowitz minimum-variance set              ####
################################################################

### Redefine the function MinimumVarianceSet from Lab 1
MinimumVarianceSet <- function(MeanRet, CovRet, Num){
  N       <- dim(CovRet)[1]
  muPF    <- seq(from=min(MeanRet)+0.001, to=max(MeanRet)-0.001, length=Num)
  A       <- cbind(rep(1,N), MeanRet)
  sdPF    <- rep(0, Num)
  weights <- matrix(0, Num, N)
  
  for(i in 1:Num){
    b           <- c(1, muPF[i])
    sol         <- solve.QP(Dmat=2*CovRet, dvec=rep(0,N), Amat=A, bvec=b, meq=2)
    sdPF[i]     <- sqrt(sol$value)
    weights[i,] <- sol$solution
  }
  
  out         <- list()
  out$weights <- weights
  out$sdPF    <- sdPF
  out$muPF    <- muPF
  return(out)
}

### MV frontier
MV <- MinimumVarianceSet(muHat, SigHat, 2000)

### Overlay the two frontiers in one plot
plot_idx <- c(which(seq_along(MV$sdPF) <= (idx_gmv <- which.min(MV$sdPF)) &
                      MV$sdPF <= max(MV$sdPF[idx_gmv:length(MV$sdPF)])),
              (idx_gmv:length(MV$sdPF))[-1])

plot(MV$sdPF[plot_idx], MV$muPF[plot_idx], type="l", lwd=2, col="#b2e061",
     xlab="Portfolio Standard Deviation",
     ylab="Portfolio Mean Return",
     main="EU vs MV Frontier (short-sales allowed)")
lines(sd_EU, mu_EU, lwd=2, lty=2, col="#bd7ebe")
legend("bottomright", legend=c("MV","EU"),
       col=c("#b2e061","#bd7ebe"), lty=c(1,2), lwd=2)

dev.copy(png, filename="EU+MV_frontier.png", width=800, height=600)
dev.off()

################################################################
### Question 3: smallest lambda vs largest lambda vs GMV     ###
################################################################

### Optimal portfolio at the smallest and largest lambda in the grid
mu_Min <- mu_EU[1]
sd_Min <- sd_EU[1]
mu_Max <- mu_EU[NLam]
sd_Max <- sd_EU[NLam]

### GMV = point with smallest sd on the MV frontier from Question 2
iGMV   <- which.min(MV$sdPF)
mu_GMV <- MV$muPF[iGMV]
sd_GMV <- MV$sdPF[iGMV]

### Compare the portfolios
CompareTable <- data.frame(
  Portfolio = c("EU (smallest lambda)", "EU (largest lambda)", "GMV"),
  Lambda    = c(round(Lambda[1], 2), round(Lambda[NLam], 2), NA),
  Mean      = c(round(mu_Min, 4), round(mu_Max, 4), round(mu_GMV, 4)),
  SD        = c(round(sd_Min, 6), round(sd_Max, 6), round(sd_GMV, 6))
)
print(CompareTable, row.names = FALSE)

################################################################
### Question 4: GMM estimation of CCAPM and FF3              ###
### dropping Size 1 portfolios (P01,...,P05)                 ###
################################################################

### Convert all net returns and factors into gross returns
Returns1 <- Returns[, 6:25]          # P06,...,P25
Ret1_G   <- Returns1/100 + 1
RM_G     <- Dataset[, 2]/100 + 1
SMB_G    <- Dataset[, 3]/100 + 1
HML_G    <- Dataset[, 4]/100 + 1
CG_G     <- Dataset[, 6]/100 + 1

### Define the moment function g for CCAPM:
### m_t(theta) = theta0 + theta1*CG_G
gCCAPM <- function(Paras, Data){    
  N   <- dim(Data)[2]     
  R   <- as.matrix(Data[,1:(N-1)])
  CG  <- Data[,N]
  m   <- Paras[1] + Paras[2]*CG
  out <- R * as.vector(m) - 1
  return(out)
}

CCAPMData <- cbind(Ret1_G, CG_G)
CCAPMout  <- gmm(gCCAPM, CCAPMData, c(0,0))
summary(CCAPMout)

### Define the moment function g for FF3:
### m_t(theta) = theta0 + theta1*RM_G + theta2*SMB_G + theta3*HML_G
gFF3 <- function(Paras,Data){
  N   <- dim(Data)[2]
  R   <- as.matrix(Data[,1:(N-3)])
  RM  <- Data[,N-2]
  SMB <- Data[,N-1]
  HML <- Data[,N]
  m   <- Paras[1] + Paras[2]*RM + Paras[3]*SMB + Paras[4]*HML
  out <- R * as.vector(m) - 1
  return(out)
}
FF3Data  <- cbind(Ret1_G,RM_G,SMB_G,HML_G)
FF3out   <- gmm(gFF3, FF3Data, c(0,0,0,0))
summary(FF3out)

################################################################
### Summary table (point estimates, SEs, J-test)             ###
################################################################

CCAPMCoef        <- coef(CCAPMout)
names(CCAPMCoef) <- c("theta0","theta1")
CCAPMSE          <- coef(summary(CCAPMout))[,"Std. Error"]
df_CCAPM         <- CCAPMout$df
J_CCAPM          <- as.numeric(specTest(CCAPMout)$test[1])
p_CCAPM          <- as.numeric(specTest(CCAPMout)$test[2])

FF3Coef        <- coef(FF3out)
names(FF3Coef) <- c("theta0","theta1","theta2","theta3")
FF3SE          <- coef(summary(FF3out))[,"Std. Error"]
df_FF3         <- FF3out$df
J_FF3          <- as.numeric(specTest(FF3out)$test[1])
p_FF3          <- as.numeric(specTest(FF3out)$test[2])

### Create summary and print
SDFSummary <- rbind(
  data.frame(
    Model     = "CCAPM",
    Parameter = c("theta0","theta1"),
    Estimate  = as.numeric(CCAPMCoef),
    SE        = as.numeric(CCAPMSE),
    J.stat    = J_CCAPM,
    df        = df_CCAPM,
    p.value   = p_CCAPM
  ),
  data.frame(
    Model     = "FF3",
    Parameter = c("theta0","theta1","theta2","theta3"),
    Estimate  = as.numeric(FF3Coef),
    SE        = as.numeric(FF3SE),
    J.stat    = J_FF3,
    df        = df_FF3,
    p.value   = p_FF3
  )
)

SDFSummaryPrint <- SDFSummary
SDFSummaryPrint$Estimate <- round(SDFSummaryPrint$Estimate, 4)
SDFSummaryPrint$SE       <- round(SDFSummaryPrint$SE, 4)
SDFSummaryPrint$J.stat   <- round(SDFSummaryPrint$J.stat, 3)
SDFSummaryPrint$p.value  <- signif(SDFSummaryPrint$p.value, 4)
print(SDFSummaryPrint)

################################################################
### Question 5: Fitted SDF series and overlay plot           ###
################################################################

### Use estimated parameters to construct the SDF implied by CCAPM and FF3
CCAPMSDF <- CCAPMCoef[1] + CCAPMCoef[2] * CG_G
FF3SDF   <- FF3Coef[1] + FF3Coef[2] * RM_G + FF3Coef[3] * SMB_G + FF3Coef[4] * HML_G

### Store the two fitted series together with the corresponding dates
SDFs <- data.frame(Date  = as.Date(Dataset$Date),
                   CCAPM = CCAPMSDF,
                   FF3   = FF3SDF)

### Produce one time-series plot with both lines overlaid
SDFsPlot <- ggplot(data=SDFs, aes(x=Date)) +
  geom_line(aes(y=CCAPM, colour="CCAPM"), linewidth=0.6) +
  geom_line(aes(y=FF3,   colour="FF3"), linewidth=0.6) +
  scale_x_date(date_breaks="5 years", date_labels="%Y") +
  scale_colour_manual(name="",
                      breaks=c("CCAPM","FF3"),
                      values=c("#bd7ebe","#b2e061")) +
  ggtitle("Estimated SDFs (drop Size 1 portfolios)") +
  theme(legend.position="bottom") +
  xlab("Year") + ylab("SDF")

SDFsPlot

ggsave(filename="SDF_CCAPM_FF3.png",
       plot=SDFsPlot, width=7, height=4.5, units="in")