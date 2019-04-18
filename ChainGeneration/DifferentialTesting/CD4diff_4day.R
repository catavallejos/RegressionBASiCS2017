# This script runs the regression model on CD4 T cells 4 days after malaria infection.
# Data was taken from Loennberg et al.

library(BASiCS)
#setwd("/nfs/research2/marioni/Nils/BASiCS/")
setwd("~/Documents/OneDrive/Projects/SingleCell/Datasets/Regression")

chains.path <- "~/Documents/OneDrive/Projects/SingleCell/BASiCS/Chains/Regression"

#### Read in data

input <- read.table("Data/Test_Data/CD4_diff.txt", sep = "\t")

#### Read in Spike-ins
ERCC.conc <- read.table("Data/Test_Data/ERCC_malaria.txt", header=TRUE, sep = "\t")

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,2]*(10^(-18)))*(6.0221417*(10^23))*9e-3
ERCC.num.final <- ERCC.num
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,1]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
			   stringsAsFactors = FALSE)


input <- input[,grepl("4_", colnames(input))]
chips <- sapply(colnames(input), function(n){unlist(strsplit(n, "\\."))[1]})

# Generate data object
Data.4day <- newBASiCS_Data(Counts = input, Tech = grepl("ERCC", rownames(input)), 
                            SpikeInfo = SpikeInput.1, BatchInfo=chips)

#### Run MCMC on these condition

MCMC.4day <- BASiCS_MCMC(Data.4day, 80000, 40, 40000, 
                         Regression = TRUE, PrintProgress=TRUE,
                         StoreChains = TRUE, StoreDir = chains.path, 
                         RunName = "CD4diff_4day_long")

#saveRDS(MCMC.4day, "Tdist/Results/Differential_testing/MCMC_CD4diff_4day.rds")
