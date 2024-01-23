
# Install LDheatmap from GitHub
## remotes::install_github("SFUStatgen/LDheatmap")


## install GAPIT from github
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT", force=TRUE)
library(GAPIT)

##install VcfR and othere dependencies
install.packages("vcfR")
install.packages("tidyverse")

library(vcfR)
library(tidyverse)

# Set working directory to correct folder, create folder GWAS in each computer, 
#and set it as working directory

getwd()
sandeep <- getwd()
setwd(sandeep)

# Read genotypic data
# Use of package vcfR to read vcf files, extract relevant information, 
# and convert it to formats usable by GAPIT

#file reading
geno <- read.vcfR("soynam_geno_sub.vcf")
geno

#extract marker information from vcf file
geno1<- extract.gt(geno)
geno1

#show 10 rows and 10 columns, just to check
geno1[1:10,1:10]


# Converting marker formats to single digits (0 and 2 for homozygotes, 1 for heterozygotes)
geno2 <- ifelse(geno1 == "0/0", 0,
                        ifelse(geno1 == "1/1", 2, 1))

#show 10 rows and 10 columns, just to check
geno2[1:10,1:10]


#transpose matrix because GAPIT needs accessions in rows, markers in columns
geno3 <- t(geno2)


#show 10 rows and 10 columns, just to check
geno3[1:10,1:10]


#conversion of file from a matrix to a data frame, to allow easy manipulation; 
#the result is creating a 
#new column with the row names
geno4 <- as.data.frame(cbind(rownames(geno3),
                                     geno3))
geno4[1:10,1:10]

#new column, number 1, named as "taxa"
colnames(geno4)[1] <- "taxa"
geno4[1:10,1:10]



# Assuming your data frame is named 'geno_numeric2'
threshold <- 0.2  # Set the threshold for missing values (20%)

# Calculate the proportion of missing values in each column
missing_values_proportion <- colMeans(is.na(geno4), na.rm = TRUE)

# Subset the data frame to include only columns with less than or equal to 20% missing values
geno5 <- geno4[, missing_values_proportion <= threshold]

# Print the cleaned data frame
print(geno5)

geno5 [1:20, 1:20]


#in this case, marker data were read as text, not numbers; next line converts 
#marker columns to numbers
geno5[,2:ncol(geno5)] <- sapply(geno5[,2:ncol(geno5)],
                                                as.numeric)
geno5[1:15,1:15]


# Install and load the imputeTS package if not already installed
# install.packages("imputeTS")
library(imputeTS)

# Replace NA values with imputed values using linear interpolation
geno6 <- na.interpolation(geno5, option = "linear")

# Alternatively, you can use other imputation methods like 'spline', 'stine', etc.
# geno_numeric2_imputed <- na.interpolation(geno_numeric2, option = "spline")

# View the imputed data
print(geno6)

geno6[1:20, 1:20]

#creation of new file with command from vcfR, which extracts relevant 
#information for each marker:
#identifier, chromosome, position
genome_map <- getFIX(geno)
genome_map 


#reorder columns as required by GAPIT: ID, chromosome, position, and rename 
#columns as required by GAPIT
genome_map1 <- as.data.frame(genome_map[,c(3,1,2)])
genome_map1

colnames(genome_map1) <- c("Name", "Chromosome", "Position")
genome_map1

head(genome_map1)




# Read phenotypic data
#command read.csv2 reads csv filed with semi colon as column separator; path 
#may have to be adjusted

pheno_data <- read.csv("soyname_pheno_all.csv")

#selects columns with target variables, not essential
pheno_data <- pheno_data[, c(1:3, 6)]

#removes rows (genotypes) with missing data in any variable
pheno_data <- pheno_data[complete.cases(pheno_data),]

#create a new object with intersection of files geno6 and pheno_data, 
#with rows in which $id includes $taxa
geno7 <- geno6[geno6$taxa %in% pheno_data$RIL,]

#same as previous command, but opposite
pheno_data_2 <- pheno_data[pheno_data$RIL %in% geno6$taxa,]

#this is not needed, but explains the next command
summary(geno7$taxa == pheno_data_2$RIL)

#match order of files; GAPIT requires same order of genotypic and phenotypic data  
pheno_data_2 <- pheno_data_2[match(geno7$taxa, pheno_data_2$RIL),]

pheno_data_2
#change name of first column, to give the same name as recommended in the GAPIT manual
colnames(pheno_data_2)[1] <- c("Taxa")

pheno_data_2
geno7[1:15, 1:15]

# Identify common column names between genotpe file geno7 qnd map file 
common_columns <- intersect(colnames(geno7), genome_map1$Name)

# Subset genome_map1_map based on common column names
genome_map2 <- genome_map1[genome_map1$Name %in% common_columns, ]

# Print the head of the resulting data frame
print(head(genome_map2))
dim(genome_map2)


# Run GAPIT

GWAS_model <- "GLM"
trait <- "Protein"

#Create folders to store GAPIT outputs without overwritting
GWAS_folder <- file.path(sandeep, GWAS_model, trait)
dir.create(GWAS_folder, showWarnings = FALSE, recursive = TRUE)
setwd(GWAS_folder)
#running GAPIT
myGAPIT <- GAPIT(
  Y = pheno_data_2[,c(1,3)],
  GD = geno7,
  GM = genome_map2,
  SNP.MAF=0.05,
  PCA.total = 3,
  model="FarmCPU"
)


warnings()

