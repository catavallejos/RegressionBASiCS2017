#!/usr/bin/env Rscript

chains.path <- "~/Documents/OneDrive/Projects/SingleCell/BASiCS/Chains/Regression"

##########################################################
#### Script to run the model on microglia cells ##########
##########################################################

# The regression and non-regression model is run on microglia cells 
# from Zeisel et al. 
# The script takes the number of GRBFs, their scale 
# parameter and the degrees of freedom as input.

#setwd("/nfs/research2/marioni/Nils/BASiCS/")
setwd("~/Documents/OneDrive/Projects/SingleCell/Datasets/Regression")

library(BASiCS)

args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1])

Var = as.numeric(args[2])

eta = as.numeric(args[3])

#### Zeisel data ####
input.Zeisel <- read.table("Data/Test_Data/microglia_Zeisel.txt", sep = "\t")

# Read in ERCCs 
ERCC.conc <- read.table("Data/cms_095046.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.Zeisel)[grepl("ERCC", rownames(input.Zeisel))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

# Generate Data object
Data.Zeisel <- newBASiCS_Data(Counts = input.Zeisel, 
                              Tech = grepl("ERCC", rownames(input.Zeisel)), 
                              SpikeInfo = SpikeInput.1)


MCMC.Zeisel <- BASiCS_MCMC(Data = Data.Zeisel, N=20000, Thin = 10, Burn = 10000, 
                           Regression = TRUE, #k=k, Var=Var, eta=eta,
                           StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_Zeisel_new_prior")


# Run non-regression model
MCMC.Zeisel.old <- BASiCS_MCMC(Data = Data.Zeisel, N=20000, Thin = 10, Burn = 10000, 
                               Regression = FALSE, Burn = 20000, PriorDelta = "log-normal",
                               StoreChains = TRUE, StoreDir = chains.path, RunName = "MCMC_Zeisel_old_new_prior")

SpikeInput.2 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput/9e-3, stringsAsFactors = FALSE)
# Generate Data object
Data.ps.2 <- newBASiCS_Data(Counts = input.ps, 
                            Tech = grepl("ERCC", rownames(input.ps)), 
                            SpikeInfo = SpikeInput.2, BatchInfo = chips)

# Comparison
Micro_OP_R <- BASiCS_LoadChain("MCMC_Zeisel", chains.path)
Micro_OP_NR <- BASiCS_LoadChain("MCMC_Zeisel_old", chains.path)


BASiCS_showFit(Micro_OP_R)
