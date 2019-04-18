#!/usr/bin/env Rscript

chains.path <- "~/Documents/OneDrive/Projects/SingleCell/BASiCS/Chains/Regression"

##############################################################
#### Script to run the model on pool-and-split RNA ###########
##############################################################

# The regression and non-regression model is run pool-and-split RNA (Grun et al). 
# The script takes the number of GRBFs, their scale 
# parameter and the degrees of freedom as input.

#setwd("/nfs/research2/marioni/Nils/BASiCS/")
setwd("~/Documents/OneDrive/Projects/SingleCell/Datasets/Regression")

# What is this?
args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1])

Var = as.numeric(args[2])

eta = as.numeric(args[3])

library(BASiCS)

#### Pool and split data ####
input.ps <- read.table("Data/Test_Data/PoolSplit.txt", sep = "\t")
# Select split condition 
input.ps <- input.ps[,grepl("RNA2i", colnames(input.ps))]
chips <- sapply(colnames(input.ps), function(n){unlist(strsplit(n, "_"))[1]})

# Read in ERCCs 
ERCC.conc <- read.table("Data/cms_095046.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.ps)[grepl("ERCC", rownames(input.ps))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

hist(log2(SpikeInput))


# Generate Data object
Data.ps <- newBASiCS_Data(Counts = input.ps, 
                          Tech = grepl("ERCC", rownames(input.ps)), 
                          SpikeInfo = SpikeInput.1, BatchInfo = chips)

# Run the regression model
MCMC.ps <- BASiCS_MCMC(Data = Data.ps, N=20000, Thin = 10, Burn = 10000, 
                       Regression = TRUE, #, k=k, Var=Var, eta=eta)
                       StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_RNA_new_prior")

#saveRDS(MCMC.ps, paste("Tdist/Results/Testing/Datasets/MCMC_RNA2i_", 
#                       k, "_", Var, "_", eta, "_reg.rds", sep=""))

# Run the non-regression model
MCMC.ps.old <- BASiCS_MCMC(Data = Data.ps, N=20000, Thin = 10, Burn = 10000, 
                           Regression = FALSE, PriorDelta = "log-normal",
                           StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_RNA_old_new_prior")

#saveRDS(MCMC.ps.old, paste("Tdist/Results/Testing/Datasets/MCMC_RNA2i_old.rds", 
#                           sep=""))


# 

SpikeInput.2 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput/9e-3, stringsAsFactors = FALSE)
# Generate Data object
Data.ps.2 <- newBASiCS_Data(Counts = input.ps, 
                          Tech = grepl("ERCC", rownames(input.ps)), 
                          SpikeInfo = SpikeInput.2, BatchInfo = chips)

MCMC.ps <- BASiCS_MCMC(Data = Data.ps.2, N=20000, Thin = 10, Burn = 10000, 
                       Regression = TRUE, 
                       StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_RNA_new_prior_oldscale")


# Run the non-regression model
MCMC.ps.old <- BASiCS_MCMC(Data = Data.ps.2, N=20000, Thin = 10, Burn = 10000, 
                           Regression = FALSE, PriorDelta = "log-normal",
                           StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_RNA_old_new_prior_oldscale")


## Comparison
# Old prior + new ERCC
RNA_OP <- BASiCS_LoadChain("MCMC_RNA", chains.path) 
# New prior + new ERCC
RNA_NP <- BASiCS_LoadChain("MCMC_RNA_new_prior", chains.path) 
# New prior + new ERCC
RNA_NP_OS <- BASiCS_LoadChain("MCMC_RNA_new_prior_oldscale", chains.path) 

p1 <- BASiCS_diagPlot(RNA_OP)
p2 <- BASiCS_diagPlot(RNA_NP)
p3 <- BASiCS_diagPlot(RNA_NP_OS)

library(cowplot)
plot_grid(p1, p2, p3, ncol = 3)

BASiCS_showFit(RNA_OP)
BASiCS_showFit(RNA_NP)
BASiCS_showFit(RNA_NP_OS)

plot(colMedians(displayChainBASiCS(RNA_OP, Param = "mu")),
     colMedians(displayChainBASiCS(RNA_OP_OS, Param = "mu")), log = "xy")

plot(colMedians(displayChainBASiCS(RNA_NP, Param = "mu")),
     colMedians(displayChainBASiCS(RNA_NP_OS, Param = "mu")), log = "xy")

