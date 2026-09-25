############################################################
### Factor models, PCA and Fama-French                   ###
### Mathematical and Empirical Finance, 2025-2026        ###
### Run from the repository folder: Rscript src/assignment2.R
############################################################

rm(list = ls())
suppressPackageStartupMessages(library(ggplot2))

##########################################################
### Load data                                          ###
##########################################################
Dataset      <- read.csv("data/DataAssignment2.csv", check.names = FALSE)
Dates        <- as.Date(Dataset[, 1])
ExcessMarket <- Dataset[, 24]
RiskFree     <- Dataset[, 25]
SMB          <- Dataset[, 26]
HML          <- Dataset[, 27]
Returns      <- as.matrix(Dataset[, 2:23])
colnames(Returns) <- sub("RET_", "", colnames(Returns))

t <- nrow(Returns)            # 1258 daily observations
N <- ncol(Returns)            # 22 stocks
BigTechNames <- colnames(Returns)[1:6]

dir.create("figures", showWarnings = FALSE)

######################################################
### Question 1: number of principal components     ###
######################################################

### PCA on the covariance matrix of returns (all returns are in %)
PCAResult       <- prcomp(Returns, center = TRUE, scale. = FALSE)
VarExplained    <- PCAResult$sdev^2 / sum(PCAResult$sdev^2)
CumVarExplained <- cumsum(VarExplained)

ScreeData <- data.frame(PC = 1:N, Individual = VarExplained,
                        Cumulative = CumVarExplained)

ScreePlot <- ggplot(ScreeData, aes(x = PC)) +
  geom_col(aes(y = Individual), fill = "grey75") +
  geom_line(aes(y = Cumulative), linewidth = 1) +
  geom_point(aes(y = Cumulative)) +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  scale_x_continuous(breaks = 1:N) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Principal component") +
  ylab("Share of total variance") +
  ggtitle("Scree plot: individual (bars) and cumulative (line) variance explained") +
  theme_minimal()
ggsave("figures/ScreePlot.png", ScreePlot, width = 7, height = 4.5, units = "in", dpi = 200)

### K = number of components with cumulative variance closest to 90%
K <- which.min(abs(CumVarExplained - 0.90))
cat("K =", K, " cumulative variance =", round(CumVarExplained[K], 4), "\n")
print(round(CumVarExplained, 4))

#############################################################
### Question 2: PCR(K) versus FF3                         ###
#############################################################

### First K principal components as factors
Factors1 <- PCAResult$x[, 1:K, drop = FALSE]

### Excess stock returns for the FF3 model
ExcessReturns <- Returns - RiskFree

### Multivariate regressions for all 22 stocks
PCROut1 <- lm(Returns ~ Factors1)
FF3Out1 <- lm(ExcessReturns ~ ExcessMarket + SMB + HML)

### Adjusted R^2, averaged over the six big-tech stocks
AdjR2_PCR <- sapply(summary(PCROut1), function(x) x$adj.r.squared)
AdjR2_FF3 <- sapply(summary(FF3Out1), function(x) x$adj.r.squared)
names(AdjR2_PCR) <- names(AdjR2_FF3) <- colnames(Returns)
AdjR2_PCR_BT <- mean(AdjR2_PCR[1:6])
AdjR2_FF3_BT <- mean(AdjR2_FF3[1:6])

### Residual correlations for the six big-tech stocks
ResidCorr_PCR_BT <- cor(PCROut1$residuals[, 1:6])
ResidCorr_FF3_BT <- cor(FF3Out1$residuals[, 1:6])
AvgAbsResidCorr_PCR_BT <- mean(abs(ResidCorr_PCR_BT[upper.tri(ResidCorr_PCR_BT)]))
AvgAbsResidCorr_FF3_BT <- mean(abs(ResidCorr_FF3_BT[upper.tri(ResidCorr_FF3_BT)]))

Q2Table <- data.frame(
  Model = c(paste0("PCR(", K, ")"), "FF3"),
  AvgAdjR2_BigTech = c(AdjR2_PCR_BT, AdjR2_FF3_BT),
  AvgAbsResidCorr_BigTech = c(AvgAbsResidCorr_PCR_BT, AvgAbsResidCorr_FF3_BT))
print(Q2Table, digits = 4)
print(round(rbind(PCR = AdjR2_PCR, FF3 = AdjR2_FF3), 3))
print(round(ResidCorr_PCR_BT, 3))
print(round(ResidCorr_FF3_BT, 3))

######################################################################
### Question 3: implied covariance matrices                        ###
######################################################################

### Diagonal idiosyncratic covariance with model-specific dof adjustment
ResCovPCR <- (t - 1) / (t - K - 1) * diag(diag(cov(PCROut1$residuals)))
ResCovFF3 <- (t - 1) / (t - 4)     * diag(diag(cov(FF3Out1$residuals)))

### Factor loadings (drop the intercepts), N x K and N x 3
BetaPCR <- t(coef(PCROut1)[-1, , drop = FALSE])
BetaFF3 <- t(coef(FF3Out1)[-1, ])

FF3Fact1 <- cbind(ExcessMarket, SMB, HML)

CovRetPCR <- BetaPCR %*% cov(Factors1) %*% t(BetaPCR) + ResCovPCR
CovRetFF3 <- BetaFF3 %*% cov(FF3Fact1) %*% t(BetaFF3) + ResCovFF3
SampleCov <- cov(Returns)

CovRetPCR_BT <- CovRetPCR[1:6, 1:6]
CovRetFF3_BT <- CovRetFF3[1:6, 1:6]
SampleCov_BT <- SampleCov[1:6, 1:6]
dimnames(CovRetPCR_BT) <- dimnames(CovRetFF3_BT) <- dimnames(SampleCov_BT)

print(round(CovRetPCR_BT, 3))
print(round(CovRetFF3_BT, 3))
print(round(SampleCov_BT, 3))

FrobPCR <- norm(CovRetPCR_BT - SampleCov_BT, type = "F")
FrobFF3 <- norm(CovRetFF3_BT - SampleCov_BT, type = "F")
FrobeniusComparison <- c(PCR = FrobPCR, FF3 = FrobFF3)
print(round(FrobeniusComparison, 4))

##########################################################
### Question 4: first three PCs versus FF3 factors     ###
##########################################################
PCFactors3 <- PCAResult$x[, 1:3]
CorrPCFF <- cor(PCFactors3, cbind(MktRF = ExcessMarket, SMB = SMB, HML = HML))
print(round(CorrPCFF, 3))

### R^2 of each FF factor on the first three PCs jointly
R2_FF_on_PC <- sapply(list(MktRF = ExcessMarket, SMB = SMB, HML = HML),
                      function(f) summary(lm(f ~ PCFactors3))$r.squared)
print(round(R2_FF_on_PC, 3))

### Loadings of the first three PCs (the sign of a PC is arbitrary)
print(round(PCAResult$rotation[, 1:3], 3))
