# Libraries.
library(dplyr)
library(kableExtra)
library(knitr)
library(cobalt)
library(MatchIt)
library(dummies)
library(stringr)
library(MASS)
library(Matrix)
library(lme4)
library(R2WinBUGS)
library(coda)
library(foreign)
library(MCMCpack)
library(mnormt)
library(corpcor)
library(splines)
library(mgcv)
library(arm)

#Add the directory of the files
#setwd("")

# Now we load in the subset week 38 and accompanying propensity score variables.
synth.dat <- readRDS("./synth_week_38.RDS")
synth.ps.vars <- readRDS("./synth_ps_week_38.RDS")


# Source the MITSS function and the dependent
# files.
source("./mitss.R")

# Need to set these variables or plug them into the MITSS function.

# Need outcomes. We examine htn_meds2, dx_htn_pec2, 
# prog_to_severe_hyp, and HYPITAT. We will exclude htn_meds2 from the set
# as we will manually estimate with MITSS to begin. 
outcomes <- c( "dx_htn_pec2", "prog_to_severe_hyp", "HYPITAT")

# Assing treatment indicator.
treatmentInd <- synth.dat$tx

# We require the propensity score.
propScore <- synth.dat$ps2

# Choose how we balance.
# Choices include: all, treatment, control.
balanceOn <- "all"

# Choose the outcome(s) of interest.
outcomeVec <- synth.dat$htn_meds2

# Here we choose the data matrix on which the analysis is run.
# We will grab the data matrix we need later on.
dataMat <- NA

# Choose response curve; choices are one or two.
Resp.Curve <- "two"

# estimandMean is the function we will use within estimandFunc.
# It can be found within mitss_func.R which is a file that holds
# all functions used within mitss.R
estimandFunc <- estimandMean

# Number of estimands being computed. This can be generalized or looped over.
# We loop.
numEstimands <- 1

# Determine if your outcome is binary or otherwise.
isOutcomeBinary <- 1

# Choose the number of imputations
nNumImpute <- 40

# Defining the breakpoints/knots on the propensity score 
nNumSub <- 6

# Determine if there are any observation we will not use.
# Had we not applied the matching prior to beginning, we would
# use our matching vector here to determine which observations will
# be used.
obsUse <- NA 

# Define the confidence intervals.
confidence <- 0.95


# Get the data matrix in form for use in MITSS.
x <- model.matrix(~., data = synth.ps.vars)[,-1]
dataMat <- x 


# Run MITSS.
out <- MITSS(treatmentInd, propScore, balanceOn = "treatment",
             outcomeVec, dataMat, Resp.Curve="two",
             estimandFunc=estimandMean, numEstimands=1,
             isOutcomeBinary=1, nNumImpute=40, nNumSub= 6,
             obsUse=NA, confidence=0.95)

# Pull the relevant information from MITSS (ATT and ATE).
table_38 <- out[[1]]
table_38 <- table_38[1:2,]


# Apply MITSS to the rest of the outcomes.
out_list <- lapply(outcomes, function(x){
  
  outcomeVec <- subset(synth.dat, select = x) %>% pull(x)
  
  if(anyNA(outcomeVec) == TRUE){
    indices_na <- which(is.na(outcomeVec))
    outcomeVec <- outcomeVec[-indices_na]
    dataMat <- dataMat[-indices_na,]
    treatmentInd <- treatmentInd[-indices_na]
    output <- paste0("There are ", length(indices_na),
                     " subjects with missing values in outcome:", x)
    print(output)
    
  }
  
  # Run MITSS.
  out <- MITSS(treatmentInd, propScore, balanceOn = "treatment",
               outcomeVec, dataMat, Resp.Curve="two",
               estimandFunc=estimandMean, numEstimands=1,
               isOutcomeBinary=1, nNumImpute=40, nNumSub= 6,
               obsUse=NA, confidence=0.95)
  
})


for(i in 1:length(out_list)){
  # Pull the relevant information from MITSS (ATT and ATE).
  table_38_extend <- out_list[[i]][["super.pop"]]
  table_38 <- rbind(table_38, table_38_extend[1:2,])
  
}

t <- table_38[,4]
n = table_38[,6]

# Calculate t statistic.
tprobs <- t.prob(t,n)
table_38_2 <- table_38[,c(1, 3, 7, 8)]
table_38_2 <- cbind(table_38_2, "Pr(T>t)" = tprobs)


# Create a table.
kable(table_38_2, booktabs = TRUE, 
      digits = 3, caption = "Week 38: Causal Estimands") %>%
  kable_styling(latex_options=c("hold_position")) %>%
  group_rows("Hypertension requiring medication adjustment", 1, 2) %>%
  group_rows("Post-delivery diagnosis of HTN or PEC", 3, 4) %>%
  group_rows("Progression to severe disease", 5, 6) %>%
  group_rows("Modified HYPITAT Composite", 7, 8)


write.csv(table_38_2, "./table_38.csv")


