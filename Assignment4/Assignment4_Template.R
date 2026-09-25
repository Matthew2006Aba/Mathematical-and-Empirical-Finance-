# ==============================================================================
# Assignment 4 - Template
# CRR and BS Models, Option Pricing and Implied Volatility
# Mathematical and Empirical Finance, 2025-2026
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# BASELINE PARAMETERS (used in Questions 1-5)
# ------------------------------------------------------------------------------

# Set the baseline parameters here


# ==============================================================================
# QUESTION 1 - CRR European option prices
# ==============================================================================

EuroCRR <- function(N) {
  
  # Define here the function returning the price of the specified option
  
}

# Compute CRR prices on the grid of N values
N  <- ...
CRR <- numeric(length(N))    # preallocate space
for (i in 1:length(N)) {
  ...                        # fill CRR
}

# Display table of results
data.frame(N = N, CRR_Price = round(CRR, 4))


# ==============================================================================
# QUESTION 2 - Black-Scholes pricing function
# ==============================================================================

BSPrice <- function( ... ) {
  
  # Define here the function returning the BS price
}

# BS price of European option with baseline parameters
BS <- BSPrice( ... )
BS


# ==============================================================================
# QUESTION 3 - Plot CRR
# ==============================================================================

plot(x    = ... ,
     y    = ... , 
     type = "b", 
     log = "x", 
     pch = 19,
     xlab = "Number of time steps N (log scale)",
     ylab = "CRR price")
abline(h = BS, col = "red", lty = 2)
legend("topright", legend = c("CRR price", "BS price"),
       col = c("black", "red"), lty = c(1, 2), pch = c(19, NA))


# ==============================================================================
# QUESTION 4 - Digital CRR option
# ==============================================================================

# Set the tree parameters
N  <- ...
dt <- ...
u  <- ...
d  <- ...
q  <- ...
  
# Create here the stock tree (same as before)
S <- matrix(NA, N + 1, N + 1)
S[1, 1] <- S0
...
  
# Digital payoff at maturity
V <- matrix(NA, N + 1, N + 1)
...
  
# Backward induction
...
  
round(V[1, 1], 4)

# ==============================================================================
# QUESTION 5 - Path dependent option
# ==============================================================================
# CAREFUL: this is a path dependent option, the backward induction doesn't work!

# Set the tree parameters
N  <- ...
dt <- ...
u  <- ...
d  <- ...
q  <- ...

# All 2^N up/down sequences
Paths <- ...

# Preallocate output vectors
Payoffs  <- numeric(2^N)  
PathProb <- numeric(2^N)

# Walk each path
for (i in 1:(2^N)) {
  
  S_path <-  ...     # Preallocate values of stock path starting from S0
  
  ...  # Fill price path
  
  S_ast      <- ...               # define S* (S_0 excluded)
  Payoffs[i] <- max( ... , 0)
  
  n_up        <- ...           # total up-moves in path i
  PathProb[i] <- ...           # probability of the path i
}

# Option price = discounted risk-neutral expected payoff
PathDepPrice <- ...
round(PathDepPrice,4)

# Compare with European option (same N = 12 steps)
EuroPrice12 <- ...
round(EuroPrice12,4)


# ==============================================================================
# QUESTION 6 - Implied Volatility
# ==============================================================================

OptionData <- read.csv("Assignment_4_dataset.csv")
n_opt      <- nrow(OptionData)         # number of options

ImpliedVol <- function( ... ) {
  
  # Define here the function
  
}

# Loop over the options
IVs <- numeric(n_opt)
for (i in 1:n_opt) {
  IVs[i] <- ...
}

OptionData$IV <- IVs

# Table of the first six options
ResultTable <- data.frame(
  Strike   = OptionData$Strike,
  Type     = OptionData$Type,
  MidPrice = OptionData$MidPrice,
  IV_pct   = round(100 * OptionData$IV, 2)
)
head(ResultTable)

# Plot IV (in %) vs strike (sorted in increasing strike order)
OptionData <- ...      # sort the dataset
plot(x = ... , 
     y = ... ,
     type = "b", 
     pch = 19,
     xlab = "Strike", 
     ylab = "Implied volatility (%)",
     main = "Implied Volatility curve")