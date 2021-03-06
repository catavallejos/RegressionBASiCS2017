#!/usr/bin/env Rscript

####################################
#### Script to downsample cells ####
####################################

# The regression and non-regression model is run on CA1 pyramidal neurons 
# from Zeisel et al. 
# The script takes the number of cells, 
# a replicate id and a random seed as input.

args = commandArgs(trailingOnly=TRUE)

nocells = as.numeric(args[1])
iter = as.numeric(args[2])
seed = as.numeric(args[3])

library(BASiCS)

# Read in data
setwd("/nfs/research2/marioni/Nils/BASiCS/")

input <- read.table("Data/Test_Data/CA1_Zeisel.txt", sep = "\t")

set.seed(seed)
input <- input[,sample(1:ncol(input), nocells)]

# Read in data and calculate absolute number od spike-in molecules
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

# Generate Data object
Data <- newBASiCS_Data(Counts = input, Tech = grepl("ERCC", rownames(input)), 
                       SpikeInfo = SpikeInput.1)

# Run the regression model
MCMC <- BASiCS_MCMC(Data = Data, 40000, 20, 20000, PrintProgress=FALSE, 
                    Regression=TRUE, PriorDelta="log-normal")

saveRDS(MCMC, paste("Tdist/Results/Testing/Downsampling_Zeisel/", "MCMC_",
                    nocells, "_", iter, ".rds", sep=""))

# Run the non-regression model
MCMC <- BASiCS_MCMC(Data = Data, 40000, 20, 20000, PrintProgress=FALSE, 
                    PriorDelta="log-normal")

saveRDS(MCMC, paste("Tdist/Results/Testing/Downsampling_Zeisel_old/", "MCMC_", 
                    nocells, "_", iter, ".rds", sep=""))



