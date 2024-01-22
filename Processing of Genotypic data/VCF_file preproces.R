install.packages("vcfR")
install.packages("tidyverse")

library(vcfR)
library(tidyverse)

# Set working directory to correct folder, create folder GWAS in each computer, 
#and set it as working directory

get
sandeep <- getwd()
setwd(sandeep)

# Read genotypic data
# Use of package vcfR to read vcf files, extract relevant information, 
# and convert it to formats usable by GAPIT

#file reading
geno <- read.vcfR("soynam_geno_sub.vcf")
geno

#extract marker information from vcf file
geno_numeric <- extract.gt(geno)
geno_numeric

#show 10 rows and 10 columns, just to check
geno_numeric[1:10,1:10]


# Converting marker formats to single digits (0 and 2 for homozygotes, 1 for heterozygotes)
geno_numeric1 <- ifelse(geno_numeric == "0/0", 0,
                        ifelse(geno_numeric == "1/1", 2, 1))

#show 10 rows and 10 columns, just to check
geno_numeric1[1:10,1:10]


#transpose matrix because GAPIT needs accessions in rows, markers in columns
geno_numeric1 <- t(geno_numeric1)





#show 10 rows and 10 columns, just to check
geno_numeric1[1:10,1:10]

#conversion of file from a matrix to a data frame, to allow easy manipulation; 
#the result is creating a 
#new column with the row names
geno_numeric2 <- as.data.frame(cbind(rownames(geno_numeric1),
                                     geno_numeric1))
geno_numeric2[1:10,1:10]

#new column, number 1, named as "taxa"
colnames(geno_numeric2)[1] <- "taxa"
geno_numeric2[1:10,1:10]

#in this case, marker data were read as text, not numbers; next line converts 
#marker columns to numbers
geno_numeric2[,2:ncol(geno_numeric2)] <- sapply(geno_numeric2[,2:ncol(geno_numeric2)],
                                                as.numeric)
geno_numeric2[1:15,1:15]



#extract accession name from the "taxa" column, from character 1 until the end
#(specific for these data)
geno_numeric2$taxa <- substr(geno_numeric2$taxa,
                             1,
                             nchar(geno_numeric2$taxa))
geno_numeric2



#creation of new file with command from vcfR, which extracts relevant 
#information for each marker:
#identifier, chromosome, position
gendibar_map <- getFIX(geno)
gendibar_map 

#reorder columns as required by GAPIT: ID, chromosome, position, and rename 
#columns as required by GAPIT
gendibar_map <- as.data.frame(gendibar_map[,c(3,1,2)])
gendibar_map 

colnames(gendibar_map) <- c("Name", "Chromosome", "Position")
gendibar_map

#fill "Name" in file gendibar_map, which was empty, with the colnames of 
#gendibar_numeric, skipping the first column (taxa)
gendibar_map$Name <- colnames(geno_numeric2)[-1]
gendibar_map

#GAPIT requires numeric codes for chromosomes; transform codes if needed
gendibar_map$Chromosome <- 2
gendibar_map



# Read phenotypic data
#command read.csv2 reads csv filed with semi colon as column separator; path 
#may have to be adjusted

pheno_data <- read.csv("soyname_pheno_all.csv")

#selects columns with target variables, not essential
pheno_data <- pheno_data[, c(1:3, 6)]

#removes rows (genotypes) with missing data in any variable
pheno_data <- pheno_data[complete.cases(pheno_data),]

#create a new object with intersection of files gendibar_numeric and pheno_data, 
#with rows in which $id includes $taxa
geno_numeric2 <- geno_numeric2[geno_numeric2$taxa %in% pheno_data$RIL,]

#same as previous command, but opposite
pheno_data_2 <- pheno_data[pheno_data$RIL %in% geno_numeric2$taxa,]

#this is not needed, but explains the next command
summary(geno_numeric2$taxa == pheno_data_2$RIL)

#match order of files; GAPIT requires same order of genotypic and phenotypic data  
pheno_data_2 <- pheno_data_2[match(geno_numeric2$taxa, pheno_data_2$RIL),]

pheno_data_2
#change name of first column, to give the same name as recommended in the GAPIT manual
colnames(pheno_data_2)[1:2] <- c("Taxa", "Row_type")

pheno_data_2
geno_numeric2[1:15, 1:15]