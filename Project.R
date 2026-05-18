# Bio-Statistics Project
# Jorge Suazo-Victoria
rm(list=ls())
 
library(snpStats)

# These are the files on MY computer. You will need to change the file paths
bedFile <- path.expand("~/Documents/Bio_Statistics/Lab1Files/GWAS_data.bed")
bimFile <- path.expand("~/Documents/Bio_Statistics/Lab1Files/GWAS_data.bim")
famFile <- path.expand("~/Documents/Bio_Statistics/Lab1Files/GWAS_data.fam")

# step 2: read in the local files using read.plink() in snpStats
data <- read.plink(bedFile,bimFile,famFile)

clinicalFile <- "~/Documents/Bio_Statistics/Clinical_Urate.csv"
clinical<- read.csv(clinicalFile, colClasses = c("character", "numeric","factor", rep("numeric", 2), "factor", rep("numeric",2)))


# Data dictionary for the PennCath data
# FamID, The unique identifier for each individual in the study
# age, self reported age at time of enrollment
# sex, self reported biological sex (1 = male, 2, = female)
# alcohol_drinks_pw, Alcohol drinks per week  
# BMI_cov, Body mass index
# diuretic_use, use of diuretics (0 = no, 1 = yes)
# phys_activity_min_pw, physical activity minutes per week
# Serum_Urate_mgdL, Amount of urate in blood (mg/dL)

# what class is "data"
class(data)
## [1] "list"

# how many elements are in "data"
length(data)
## [1] 3

rownames(clinical) <- clinical$FamID
print(head(clinical))

clinical <- clinical[
  !(is.na(clinical$alcohol_drinks_pw) & 
      is.na(clinical$BMI_cov) & 
      is.na(clinical$diuretic_use) & 
      is.na(clinical$phys_activity_min_pw) &
      is.na(clinical$Serum_Urate_mgdL)),
]

clinical <- subset(clinical, Serum_Urate_mgdL > 0)

for (i in colnames(clinical)) {
  print(i)
  print(summary(clinical[[i]]))
}

genotype <- data$genotypes
genoMap <- data$map
colnames(genoMap)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")

snpsum.col <- col.summary(genotype)
head(snpsum.col)


## Call rate

# To filter by call rate, we remove individuals that are missing genotype data across
# more than a pre-defined percentage of the typed SNPs. This proportion of missingness
# across SNPs is the sample call rate and we apply a threshold of 95%, meaning that we
# drop those individuals missing genotype data for more than 95% of the typed SNPs. A
# new reduced dimesion SnpMatrix gentotype is created from this filter.

# Recalculate SNP summary statistics
snpsum.col <- col.summary(genotype)

# SNP call rate threshold
call.thresh <- 0.95

# Identify variants passing call rate QC
keep.call <- with(
  snpsum.col,
  !is.na(Call.rate) & Call.rate >= call.thresh
)

keep.call[is.na(keep.call)] <- FALSE

# Extract low-quality SNPs
low_call_snps <- rownames(snpsum.col)[!keep.call]

cat(
  length(low_call_snps),
  "SNPs flagged for low call rate (<95%)\n"
)

# Minor allele frequency
# Inadequate power to infer a statiscally significant relationship between
# the SNP and the trait under study is the result of a large degree of
# homogeneity at a given SNP across study participants. This occurs when 
# we have a very small minor allele frequency (MAF), which means the majority
# of the individuals have two copies of the same major alleles. Minor allele
# frequency is defined as frequency of the less common allele at a variable 
# site. In this module, we remove SNPs whose minor allele fequency is less
# than 1%. In other cases, particularly when sample size is small, we may
# apply a cut point of 5%.

minor <- 0.01

keep.minor <- with(snpsum.col, (!is.na(MAF) & MAF > minor) )
keep.minor[is.na(keep.minor)] <- FALSE


library(ggplot2)
library(dplyr)

minor <- 0.01

snpsum.col <- snpsum.col %>%
  mutate(MAF_status = ifelse(is.na(MAF), NA,
                             ifelse(MAF > minor, "Retained", "Removed")))

maf_plot <- ggplot(snpsum.col, aes(x = MAF, fill = MAF_status)) +
  geom_histogram(
    bins = 100,
    color = "white",
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = minor,
    linetype = "dashed",
    linewidth = 1.2,
    color = "red"
  ) +
  annotate(
    "text",
    x = minor,
    y = 20000,
    label = "MAF = 0.01",
    vjust = 1.5,
    color = "red",
    angle = 90
  ) +
  scale_fill_manual(
    values = c(
      "Retained" = "#3a86ff",
      "Removed" = "#d90429"
    ),
    na.value = "grey80"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top") +
  labs(
    title = "Distribution of Minor Allele Frequency",
    x = "MAF",
    y = "Number of SNPs",
    fill = ""
  )

print(maf_plot)

cat(ncol(genotype)-sum(keep.minor),"SNPs will be removed due to low MAF (",ncol(genotype),"-",sum(keep.minor),").\n"  )

# After knowing the SNPs needed to be removed, we subset the genotype and
# SNP summary data for the SNPs that pass MAF criteria.

genotype <- genotype[,keep.minor] 
snpsum.col <- snpsum.col[keep.minor,] 


# Saving work for following labs…
# Now that all the data is read into R, we need to save the objects for the
# rest of the labs running the following command:

save(genotype,genoMap,snpsum.col,genoMap,bedFile,bimFile,famFile, clinical, maf_plot, low_call_snps, file= "project1_save.RData")

#### PART 2 ====================================================================

rm(list=ls()) #remove any data objects or functions to avoid collisions or misspecification

library(snpStats)
library(SNPRelate)
library(gdsfmt)
library(dplyr)
library(ggplot2)

setwd("/home/Suaria/Documents/Bio_Statistics")
load("project1_save.RData")

#specificy the number of threads for IBD and PCA functions
numCores=4

# calculate summary statistics for the genotype matrix per variant (column)
# and per individual (row)
snpsum.row <- row.summary(genotype) #from snpStats
snpsum.col<-col.summary(genotype)


MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)

## Heterozygosity

# Heterozygosity refers to the presence of each of two alleles at a given SNP
# within an individual. This is expected under Hardy-Weinberg Equilibrium (HWE)
# to occur with probability of 2p(1-p) where p is the dominant allele frequency
# at that SNP (assuming a bi-allelic SNP). Excess or deficient heterozygosity
# across typed SNPs within an individual is commonly considered as an indication
# of inbreeding but also serves as a screen for poor sample quality.


# In this lab, heterozygosity F statistic calculation is carried out with the form,
# |F|=(1-O/E), where O is observed proportion of heterozygous genotypes for a given
# sample and E is the expected proportion of heterozygous genotypes for a given
# sample based on the minor allele frequency across all non-missing SNPs for a given
# sample. We will remove any sample with inbreeding coefficient |F| > 0.05, which
# might indicate autozygosity (which is not inherently a problem, but perhaps
# unusual) or poor quality sample.

# What type of poor quality sample is heterozygosity based quality control
# likely to identify? Cross-contamination Samples

hetExp <- callmatrix %*% (2*MAF*(1-MAF)) #Note that %*% is the 'multiply_by_matrix' function from plyr
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

plot(hetObs,hetExp);abline(0,1)
hist(snpsum.row$hetF,main="Distribution of observed vs. expected heterozygosity across variants")
head(snpsum.row)

het.thresh <- 0.05

library(ggplot2)

Distribution_of_Heterozygosity <- ggplot(snpsum.row, aes(x = hetF)) +
  
  geom_histogram(
    bins = 50,
    fill = "#4361ee",
    color = "white",
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = c(-het.thresh, het.thresh),
    color = "#d00000",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = -het.thresh,
    y = Inf,
    label = "-0.05",
    vjust = 2,
    color = "#d00000"
  ) +
  
  annotate(
    "text",
    x = het.thresh,
    y = Inf,
    label = "0.05",
    vjust = 2,
    color = "#d00000"
  ) +
  
  theme_minimal(base_size = 16) +
  
  labs(
    title = "Distribution of Heterozygosity F-statistics",
    x = "F-statistic",
    y = "Number of Individuals"
  )

keep.het <- with(snpsum.row, abs(hetF) <= het.thresh)
keep.het[is.na(keep.het)] <- FALSE  
cat(nrow(genotype)-sum(keep.het), "individuals will be removed due to low sample call rate or unusual F coefficient.\n")




# As before, after identifying the individuals to be removed, we subset the genotype
# and clinical data set using information obtained above to remove the individuals
# that do not pass the criteria from our sample.

genotype <- genotype[keep.het,]
sum(!rownames(genotype) %in% clinical$FamID)
clinical <- clinical[match(rownames(genotype), clinical$FamID), ]
# After removing individuals, summary statistics of the genetic data change.
# Recalculate sumstats
snpsum.row<-row.summary(genotype)


## Call rate

# To filter by call rate, we remove individuals that are missing genotype data across
# more than a pre-defined percentage of the typed SNPs. This proportion of missingness
# across SNPs is the sample call rate and we apply a threshold of 50%, meaning that we
# drop those individuals missing genotype data for more than 50% of the typed SNPs. A
# new reduced dimesion SnpMatrix gentotype is created from this filter.

samplecall.thresh <- 0.97

# Subset the genotype and clinical data yet again
keep.samplecall <- with(snpsum.row, !is.na(Call.rate) & Call.rate > samplecall.thresh)
keep.samplecall[is.na(keep.samplecall)] <- FALSE 
cat(nrow(genotype)-sum(keep.samplecall), "individuals will be removed due to low sample call rate.\n")

library(ggplot2)
library(dplyr)


snpsum.row <- snpsum.row %>%
  mutate(
    CallRate_Status = ifelse(
      Call.rate > samplecall.thresh,
      "Retained",
      "Removed"
    )
  )

callrate_plot <- ggplot(
  snpsum.row,
  aes(x = Call.rate, fill = CallRate_Status)
) +
  
  geom_histogram(
    bins = 50,
    color = "white",
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = samplecall.thresh,
    linetype = "dashed",
    linewidth = 1.2,
    color = "red"
  ) +
  
  annotate(
    "text",
    x = samplecall.thresh,
    y = 30,
    label = "97%",
    angle = 90,
    vjust = 1.5,
    color = "red"
  ) +
  
  scale_fill_manual(
    values = c(
      "Retained" = "#3a86ff",
      "Removed" = "#d90429"
    )
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    legend.position = "top"
  ) +
  
  labs(
    title = "Distribution of Sample Call Rate",
    x = "Sample Call Rate",
    y = "Number of Individuals",
    fill = ""
  )

print(callrate_plot)

genotype <- genotype[keep.samplecall,]
clinical <- clinical[match(rownames(genotype), clinical$FamID), ]
snpsum.row<-row.summary(genotype)


# Hardy-Weinberg Equilibrium (HWE)
# In case-control studies, HWE is typically evaluated in controls only,
# since disease-associated variants may deviate from equilibrium in cases.
# Because serum urate is a continuous trait and no true control group exists,
# I define a pseudo-control subset using individuals within the central
# distribution (1st-99th percentile) of serum urate values.

# Why calculate HWE based only on the CAD controls? because the affected samples don't follow the equlibrium and the filter will wipe the SNPs that we want to analyse
# How many variants are omitted comparing controls only vs. the whole sample? 0

q <- quantile(
  clinical$Serum_Urate_mgdL,
  probs = c(0.01, 0.99),
  na.rm = TRUE
)

controls <- subset(
  clinical,
  Serum_Urate_mgdL >= q[1] &
    Serum_Urate_mgdL <= q[2]
)
genotype_controls <- genotype[rownames(genotype) %in% controls$FamID, ]

n_total <- ncol(genotype)
n_controls <- ncol(genotype_controls)

var_total <- apply(genotype, 2, var, na.rm = TRUE)
var_controls <- apply(genotype_controls, 2, var, na.rm = TRUE)

omitted <- sum(var_controls == 0 & var_total != 0, na.rm = TRUE)

hardy.thresh <- 10^-6     #the threshold is parameterized as a p-value, not a z score as in col.summary
controls <- as.character(
  clinical[
    clinical$Serum_Urate_mgdL >= q[1] &
      clinical$Serum_Urate_mgdL <= q[2],
    "FamID"
  ]
) 



# use col.summary on genetic data from controls only
genotype_controls <- genotype[
  rownames(genotype) %in% controls,
]
snpsum.colControls <- col.summary(genotype_controls)
#use qnorm to convert the p value threshold 'hardy' to a z score
keep.hwe <- with(snpsum.colControls, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy.thresh/2) ) ) )

hwe_p <- 2 * pnorm(-abs(snpsum.colControls$z.HWE))

hwe_df <- data.frame(HWE_p = hwe_p)

hwe_df <- hwe_df %>%
  filter(
    !is.na(HWE_p),
    HWE_p > 0
  )

cuartiles <- quantile(hwe_df$HWE_p, probs = c(0.25, 0.5, 0.75))
tabla_cuantiles <- table(cut(hwe_df$HWE_p, 
                             breaks = c(-Inf, cuartiles, Inf), 
                             labels = c("Q1", "Q2", "Q3", "Q4")))



hwe_plot <- ggplot(
  hwe_df,
  aes(x = HWE_p)
) +
  
  geom_histogram(
    bins = 50,
    fill = "#d90429",
    color = "white"
  ) +
  
  geom_vline(
    xintercept = hardy.thresh,
    color = "blue",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  

  theme_minimal(base_size = 16) +
  
  labs(
    title = "Distribution of HWE p-values",
    x = "HWE p-value (log scale)",
    y = "Number of SNPs"
  )
# We can delete this object now
rm(snpsum.colControls)

keep.hwe[is.na(keep.hwe)] <- FALSE          

cat(ncol(genotype)-sum(keep.hwe),"SNPs will be removed due to deviation from HWE.\n") 

# Subset the genotype data yet again
genotype <- genotype[,keep.hwe]

# Check the dimensionality of the working genotype matrix
print(genotype)

## Cryptic relatedness, duplicates and gender identity

# Population-based cohort studies often only concerns unrelated individuals and
# the generalized linear modeling approach assumes independence across individuals.
# In regional cohort studies (e.g. hospital-based cohort studies) of complex
# diseases, however, individuals from the same family could be unintentionally
# recruited. To mitigate this problem, we thus employ Identity-by-descent (IBD)
# analysis which is a common measurement of relateness and duplication. Note that
# gender identity can also be checked at this stage to confirm self-reported gender
# consistency on X and Y chromosomes. With this data, we do not check gender
# identity due to unavailibility of sex chromosomes information.

# To filter on relateness criteria, we use the SNPRelate package to perform
# identity-by-descent (IBD) analysis. This package requires that the data be
# transformed into a GDS format file. IBD analysis is performed on only a subset
# of SNPs that are in linkage equilibrium by iteratively removing adjacent SNPs
# that exceed an linkage disequilibrium (LD) threshold in a sliding window using
# the snpgdsLDpruning() function. 

# Variants with pairwise linkage disequilibrium exceeding r^2 = 0.2
# were iteratively removed prior to IBD estimation.

ld.thresh <- 0.2   #maximum r2 between variants in 0.2
kin.thresh <- 0.04 #for context, the kinship coefficient for 1st cousins is 0.0625
# twins is 1, unrelated near 0

# Create a specialized 'genetic data set' or .gds object using a function from
# SNPRelate package. It reads from the original files, and subsets to the individuals and
# and variants in the current 'genotype' object
snpgdsBED2GDS(bedFile, famFile, bimFile,"GWAS_data.gds") 
genofile <- snpgdsOpen("GWAS_data.gds", readonly = FALSE)
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)


# Generate a subset of variants in low linkage disequilibrium for IBD and PCA,
# as these methods assume independence of variants
set.seed(1000)
snpSUB <- snpgdsLDpruning(genofile,
                          verbose=TRUE,
                          ld.threshold = ld.thresh,
                          maf=0.05,
                          missing.rate=0.05,
                          sample.id = rownames(genotype), 
                          snp.id = colnames(genotype)) 

snpset.LoLD <- unlist(snpSUB, use.names=FALSE)


cat(length(snpset.LoLD),"will be used in IBD analysis\n")

# Use another function in SNPRelate to estimate Identity by Descent between
# all participants using a Method of Moments estimator
geno.sample.ids <- rownames(genotype)
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE, sample.id = geno.sample.ids, snp.id = snpset.LoLD,  num.thread = numCores)
ibdcoeff_all <- snpgdsIBDSelection(ibd)    

head(ibdcoeff_all)

#select those pairs that have kinship greater than the threshold value for
# further evaluation
kins_all <- ibdcoeff_all
kins<-subset(ibdcoeff_all,ibdcoeff_all$kinship>kin.thresh) 
#ibdcoeff <- ibdcoeff_all[ ibdcoeff_all$kinship >= kin.thresh, ]

kins_all <- kins_all %>%
  mutate(
    relationship = case_when(
      kinship > 0.354 ~ "Duplicates / Twins",
      kinship > 0.177 ~ "1st Degree",
      kinship > 0.0884 ~ "2nd Degree",
      kinship > kin.thresh ~ "Related",
      TRUE ~ "Unrelated"
    )
  )

sum(duplicated(t(apply(kins[,c("ID1","ID2")], 1, sort)))) # If greater tha 0 there are duplicates

related.samples <- NULL

library(ggplot2)

ibd_plot <- ggplot(
  kins_all,
  aes(
    x = k0,
    y = k1,
    color = relationship
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  theme_minimal(base_size = 16) +
  
  scale_color_manual(
    values = c(
      "Duplicates / Twins" = "#d00000",
      "1st Degree" = "#e85d04",
      "2nd Degree" = "#ffba08",
      "Related" = "#3a86ff"
    )
  ) +
  
  labs(
    title = "Identity-by-Descent Relatedness",
    x = "IBD0 (k0)",
    y = "IBD1 (k1)",
    color = "Relationship"
  ) +
  
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0,1)
  )

#remove one member of each related pair from the dataset

related.samples <- c()
for (i in 1:nrow(kins)) {
  # Select ID based on highest k1 value
  selected_id <- ifelse(runif(1)>0.5,kins$ID1[i], kins$ID2[i])
  related.samples <- c(related.samples, selected_id)
}


#Subset the genotype and clinical data yet again
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

cat(length(related.samples), " samples removed due to kinship >=", kin.thresh,"\n") 


#evaluate the distribution of pairwise kinships
#this is ALL PAIRS of individuals, mostly minimally related (unrelated)
#plot(ibdcoeff_all$k0,ibdcoeff_all$k1,xlim=c(0,1),ylim=c(0,1),xlab="IBD0",ylab="IBD1",main="IBD of all individuals")
plot(kins$k0,kins$k1,xlim=c(0,1),ylim=c(0,1),xlab="IBD0",ylab="IBD1",main="IBD of related individuals")


## PCA

# Principle components of Ancestry is one approach to visualizing and classifying
# individuals into ancestry groups based on their observed genetic similarity.
# There are two reasons for doing this. First, self-reported ancestry can
# differ from clusters of individuals that are based on genetics infomation.
# Second, the presence of an individual not appearing to fall within a genetic
# similarity cluster may be suggestive of a sample error.

# Use another function from SNPRelate to generate PCs, by default, this function estimates
# the first 32 PCs

geno.sample.ids <- rownames(genotype)

clinical <- clinical[!is.na(clinical$FamID), ]

pca <- snpgdsPCA(
  genofile,
  sample.id = geno.sample.ids,
  snp.id = snpset.LoLD,
  num.thread = numCores,
  eigen.cnt = 0 # eigen.cnt	
  # output the number of eigenvectors; if eigen.cnt <= 0, then return all eigenvectors
)

pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10], stringsAsFactors = FALSE)

colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    
                    PC2 = pca$eigenvect[,2],
                    PC3 = pca$eigenvect[,3],
                    PC4 = pca$eigenvect[,4],
                    PC5 = pca$eigenvect[,5],
                    PC6 = pca$eigenvect[,6],
                    stringsAsFactors = FALSE)

plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot",pch="o")

pctab$Urate <- clinical$Serum_Urate_mgdL[
  match(pctab$sample.id, clinical$FamID)
]


pctab <- pctab[!is.na(pctab$Urate), ]

library(plotly)

pca_3d <- plot_ly(
  pctab,
  
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  
  color = ~Urate,
  
  type = "scatter3d",
  mode = "markers",
  
  marker = list(
    size = 3,
    opacity = 0.8
  )
) %>%
  
  layout(
    title = "3D PCA of Genetic Structure",
    
    scene = list(
      xaxis = list(
        title = paste0(
          "PC1 (",
          round(pca$varprop[1] * 100, 2),
          "%)"
        )
      ),
      
      yaxis = list(
        title = paste0(
          "PC2 (",
          round(pca$varprop[2] * 100, 2),
          "%)"
        )
      ),
      
      zaxis = list(
        title = paste0(
          "PC3 (",
          round(pca$varprop[3] * 100, 2),
          "%)"
        )
      )
    )
  )

pca_3d

pca_plot <- ggplot(pctab, aes(PC1, PC2 , color=Urate)) +
  geom_point() +
  stat_ellipse()

# What about the other PCs? Any structure apparent there?
# To evaluate the pattern of how much genetic variation is captured by PCs,
# we simply plot the eigenvalue.

eigen <- plot(pca$eigenval,ylab="Eigenvalue",xlab="Index of Eigenvector",xlim=c(1,20),pch=20,main="Scree plot")

df <- data.frame(
  PC = seq(length(pca$varprop)),
  variance = pca$eigenval
)

df <- df[1:10,]

pvar_plot <- ggplot(df, aes(x = PC, y = variance)) +
  geom_col(fill = "#ff7733") +           # Bar plot for individual variance
  geom_line() +                            # Line to show the "elbow"
  geom_point(size = 3) + 
  scale_x_continuous(breaks = 1:nrow(df)) +
  labs(title = "Scree Plot",
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  theme_minimal()

# What would we need to change to estimate all PCs from which we could find
# the trace and express the eigenvalues as the proportion of variance?
# the snpgdsPCA captures this measure for the the calculated eigenvalues, check out pca$varprop 

plot(pca$varprop,ylab="Total Varianza",xlab="Index of Eigenvector",xlim=c(1,20),pch=20,main="Scree plot")


# Terminate our special .gds file. R gets cranky if you close a session and
# the link to the gds is open
closefn.gds(genofile)

#save the processing in this step for the next lab
#write.csv(pctab, file = "pctab.csv", row.names = FALSE)
save(genotype,genoMap,snpsum.col, snpsum.row, genofile, clinical, pcs, Distribution_of_Heterozygosity, hwe_plot, callrate_plot, ibd_plot, pca_plot, pvar_plot, pca_3d, low_call_snps, file= "project2_save.RData")

# PART 3 =======================================================================

rm(list=ls()) #remove any data objects or functions to avoid collisions or misspecification

#BiocManager::install("doParallel",force=TRUE) #run only once to install
#install.packages("devtools") #run only once to install

# module load HGI/softpack/groups/rnamelio/pb20-R_wgs_bulkRNA/4
# 

library(dplyr)
library(doParallel)
library(devtools)

source("https://raw.githubusercontent.com/CrisVanHout-LIIGH/Class/refs/heads/main/LabFunctions.R")

setwd("~/Documents/Bio_Statistics")
load("project2_save.RData")  

#Print the number of SNPs to be checked
cat(paste(ncol(genotype), "Thus far, SNPs are included in analysis.\n"))
cat(paste(nrow(genotype), "Thus far, participants are included in analysis.\n"))

if(file.exists("GWAA_final.txt")){
  
  cat("Previous GWAS detected...\n")
  
  # Cargar el archivo final
  gwas <- read.table("GWAA_final.txt", header = FALSE)
  colnames(gwas) <- c("SNP", "Estimate", "Std.Error", "t.value", "p.value")
  
  cat(nrow(gwas), "Total Processed SNPs\n")
  
  # skip
  skip_gwas <- TRUE
  
} else {
  
  cat("Starting GWAS...\n")
  skip_gwas <- FALSE
  
}

#Manipulating Phenotype Data
# Next we focus on the phenodata argument. First we assign phenoSub to be 'phenodata',
# which is created from merging the clinical and pcs data. We then create the phenotype
# value (dependent or outcome variable) by choosing (or not) a transformation of the phenotype
# variable. Next we remove anyone that doesn’t have data available, and subset the genotype
# data for the remaining individuals. Lastly, we want to remove variables that won’t be used
# in analysis, e.g. the remaining outcome variables and the untransformed version of the phenotype.

traitName <- "Serum_Urate_mgdL"
hist(clinical$Serum_Urate_mgdL)

clinicalSub <- clinical[, c(
  "FamID",
  "sex",
  "age",
  "alcohol_drinks_pw",
  "BMI_cov",
  "diuretic_use",
  "phys_activity_min_pw",
  traitName
)]

names(clinicalSub)[names(clinicalSub) == traitName] <- "raw_trait"
clinicalSub$phenotype <- clinicalSub$raw_trait # or log(), rntransform(), apply a transformation here
clinicalSub$raw_trait <- NULL


# Subset pcs to keep FamID and any PCs that you choose to include
pcsSub <- pcs[, c("FamID", "pc1", "pc2", "pc3", "pc4")]


# Merge clinicalSub and pcsSub by FamID
phenodata <- merge(clinicalSub, pcsSub, by = "FamID")

# Remove rows with missing phenotype values
phenodata <- phenodata[!is.na(phenodata$phenotype), ]

# Align genotype rows to phenodata
phenodata <- phenodata[phenodata$FamID %in% rownames(genotype), ]
genotype  <- genotype[phenodata$FamID, ]

epidemiological_model <- lm(
  phenotype ~ . - FamID +
    BMI_cov*phys_activity_min_pw,
  data = phenodata
) #not testing any snps, but all covariates included in phenodata (except FamID)
summary(epidemiological_model)


# Note that the GWAS model in the GWAA function takes all variables in the 'phenodata' data.frame as
# covariates, i.e. age sex and any pcs. One may control which covariates are included in
# the phenodata data.frame, or by changing the glm in the GWAA function to explicitly account
# for covariates


#The GWAA function
#The function itself has four major elements. Given the data described above, it
# (1) determines how many cores to use/which method to employ for parallel processing
# (2) converts the genotype data into the necessary format for an additive model,
# (3) fits a linear model for each SNP and covariates (one group at a time) and then
# (4) writes the SNP coefficient details (the estimate, standard error, test statistic
# and p value) from each model out to the GWAA.txt file.

# Note that each core uses about 250Mb memory to analyze all autosomes ~350k variants, or 2Gb total.
# Low memory computers may want to reduce the number of hosts to 2 or even 1 to limit use of RAM.
# (and/or close your memory hungry web browser while you run this analysis!)

GWAA <- function(genotype, phenodata, filename = NULL,
                 append = FALSE,
                 workers = getOption("mc.cores", 2L),
                 flip = FALSE,
                 select.snps = NULL,
                 hosts = 4) {
  
  if (is.null(hosts)) {
    cl <- makeCluster(workers)
  } else {
    cl <- makeCluster(hosts, "PSOCK")
  }
  
  show(cl)
  registerDoParallel(cl)
  
  flip.matrix <- function(x) {
    zero2 <- which(x == 0)
    two0 <- which(x == 2)
    
    x[zero2] <- 2
    x[two0] <- 0
    
    return(x)
  }
  
  foreach(part = 1:nSplits) %do% {
    
    outfile <- paste0("chunk_", part, ".txt")
    
    # Skip chunks
    if(file.exists(outfile)){
      cat("Chunk", part, "exists, skipping\n")
      NULL 
    } else {
      # TODO 
      genoNum <- as(
        genotype[, snp.start[part]:snp.stop[part]],
        "numeric"
      )
      
      if (isTRUE(flip))
        genoNum <- flip.matrix(genoNum)
      
      rsVec <- colnames(genoNum)
      
      res <- foreach(
        snp.name = rsVec,
        .combine = "rbind"
      ) %dopar% {
        
        pheno_tmp <- phenodata
        
        pheno_tmp$snp <- genoNum[, snp.name]
        
        fit <- glm(
          phenotype ~ . - FamID +
            BMI_cov:phys_activity_min_pw,
          family = gaussian,
          data = pheno_tmp
        )
        
        coef_out <- summary(fit)$coefficients["snp", ]
        
        data.frame(
          SNP = snp.name,
          Estimate = coef_out["Estimate"],
          Std.Error = coef_out["Std. Error"],
          t.value = coef_out["t value"],
          p.value = coef_out["Pr(>|t|)"]
        )
        
      }
      
      write.table(
        res,
        outfile,
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE
      )
      
      cat(sprintf(
        "GWAS SNPs %s-%s (%s%% finished)\n",
        snp.start[part],
        snp.stop[part],
        round(100 * part/nSplits, 2)
      ))
    }
  }
  # UNIR TODOS LOS CHUNKS
  files <- list.files(
    pattern = "^chunk_[0-9]+\\.txt$"
  )
  
  files <- sort(
  list.files(pattern="^chunk_[0-9]+\\.txt$")
)
  
  final <- do.call(
    rbind,
    lapply(files, read.table, header = FALSE)
  )
  
  
  
  write.table(
    final,
    "GWAA_final.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  stopCluster(cl)
  
  return("Done.")
}
#this takes about 30seconds on 8 cores, chr10 is close to the genomewide average
#SNPsToTest <- genoMap$SNP[genoMap$chr==21]
#so, this should take about 22 times as long as chr10
SNPsToTest <- genoMap$SNP[genoMap$chr<=22] 

genotype_subset <- genotype[,colnames(genotype)%in%SNPsToTest]

nSNPs <- ncol(genotype_subset)
nSplits <- 20
genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

#Run the GWAS
if(!skip_gwas){
  
  SNPsToTest <- genoMap$SNP[genoMap$chr<=22] 
  genotype_subset <- genotype[,colnames(genotype)%in%SNPsToTest]
  
  nSNPs <- ncol(genotype_subset)
  nSplits <- 20
  genosplit <- ceiling(nSNPs/nSplits)
  snp.start <- seq(1, nSNPs, genosplit)
  snp.stop <- pmin(snp.start+genosplit-1, nSNPs)
  
  start <- Sys.time()
  GWAA(
    genotype_subset,
    phenodata,
    filename="GWAA.txt",
    hosts = 5
  )
  end <- Sys.time()
  print(end-start)
  
  gwas <- read.table("GWAA_final.txt", header = FALSE)
  colnames(gwas) <- c("SNP", "Estimate", "Std.Error", "t.value", "p.value")
  
} else {
  cat("Using GWAS previous results\n")
}

gwas <- gwas[!(gwas$SNP %in% low_call_snps), ]

cat(
  nrow(gwas),
  "GWAS SNPs retained after post-hoc SNP call rate filtering\n"
)

gwasAnnotation<- gwas %>% left_join(genoMap, by = "SNP") #Annotate variant information for any SNP that exists in the gwas output. left_join annotates any SNP in the gwas data.frame with information from genoMap. right_join annotates any SNP in the genoMap data.frame with information from gwas. #note that base R 'merge' threw an error that the vectors in the merge were too long, which they shouldn't be, but this works as expected
colnames(gwasAnnotation)[colnames(gwasAnnotation) == "chr"] <- "CHR"
colnames(gwasAnnotation)[colnames(gwasAnnotation) == "position"] <- "BP"
colnames(gwasAnnotation)[colnames(gwasAnnotation) == "p.value"] <- "P"


# now, some essential visualizations
par(mfrow = c(1, 1)) 
bonCorrection<-0.05/nrow(gwas)
bonLine <- -log10(bonCorrection)
manhattan(
  gwasAnnotation,
  suggestiveline = FALSE,
  genomewideline = bonLine
)
QQ_plot(gwasAnnotation$P)

gwasAnnotation[order(gwasAnnotation$P), ][1:5, ]

# Extract a significantly associated variant from genotype matrix
variant <- gwasAnnotation[order(gwasAnnotation$P), ][1,1 ] #select top association
#variant <- "rs1532625" #or select a specific variant by rsID
snp_subset <- genotype[, c(variant), drop = FALSE]  # Keep as dataframe

#evaluate relevant annotation for the variant and recode to the original allelic coding
variantAnnotation<-gwasAnnotation[gwasAnnotation$SNP == variant, , drop = FALSE]
variantAnnotation
snpsum.col[variant, , drop = FALSE]


write.csv(gwasAnnotation,
          "gwasAnnotation.csv",
          row.names = FALSE)

# Finally, for continuous traits, we build a small dataframe of clinical and variant data for plotting
# We may select the untransformed variable for plotting
#clinical_subset <- clinical[, c(traitName), drop = FALSE]  # Ensure it remains a dataframe
# or
# Select the transformed variable that was used in the GWAS
# Do not run for binary outcomes

if(1==1){ #set to 1==0 for binary outcomes to skip this step
  clinical_subset <- clinicalSub[, "phenotype", drop = FALSE];traitName<-"phenotype"
  rownames(clinical_subset) <- clinicalSub$FamID
  
  # Merge both pheno and geno datasets by rownames (assuming sample IDs are rownames)
  # Clean up the format
  genopheno.df <- merge(clinical_subset, snp_subset, by = "row.names")  
  genopheno.df$Row.names <- NULL
  genopheno.df[[traitName]] <- as.numeric(genopheno.df[[traitName]])
  genopheno.df[[variant]] <- as.factor(as.numeric(genopheno.df[[variant]]))
  
  # Note that this recoding is only valid if the GWAA flip=FALSE,
  # if flip=TRUE, the genotype coded by 3 is always the minor allele homozygote
  genopheno.df$variant_recoded <- ifelse(genopheno.df[[variant]] == 0, NA,
                                         ifelse(genopheno.df[[variant]] == 1, 
                                                paste(variantAnnotation$A1, variantAnnotation$A1, sep = ""),
                                                ifelse(genopheno.df[[variant]] == 2, 
                                                       paste(variantAnnotation$A1, variantAnnotation$A2, sep = ""),
                                                       paste(variantAnnotation$A2, variantAnnotation$A2, sep = ""))))
  genopheno.df$variant_recoded<-as.factor(genopheno.df$variant_recoded)
  
  # Orient the model in the context of the genotype
  head(genopheno.df)
  table(genopheno.df[[variant]])
  table(genopheno.df$variant_recoded)
  
  # Filter out rows where variant is missing
  genopheno.df.nomiss <- genopheno.df[genopheno.df[[variant]] != 0, ]
  
  library(rlang) #this package allows us to interpret the variables 'variant' and 'traitName' as below
  
  # Recall, this plot uses the phenotype values as used in the GWAS.
  mixedPlot <- ggplot(genopheno.df.nomiss, aes(x = !!sym(variant), y = !!sym(traitName), fill = !!sym(variant))) +
    geom_flat_violin(alpha = 0.5, position = position_nudge(x = .2, y = 0)) +
    geom_jitter(alpha = 0.5, width = 0.15) +
    geom_boxplot(alpha = 0.5, width = .25, outlier.shape = NA) +
    theme_minimal() +
    theme(legend.position = "none")
  
  mixedPlot
}


# You have just run a Genome Wide Association Study!
# What are the units and direction of effect for the association? The lead variant showed a positive effect size (β ≈ 0.195), indicating that the T allele is associated with increased serum urate levels under an additive genetic model.
# What gene is this variant in or near? SCL2A9
# Where is that gene expressed? Chromosome 4
# What is known about the function of that gene in pathways related to the trait? Transport of urate
# Might there be diagnostic, biomarker, pharmaceutical value to this observation? There are clinical trails using SCL2A9 inhibitors

#save the processing in this step for the next lab
save(gwas,gwasAnnotation,genopheno.df,genopheno.df.nomiss,genotype,genoMap,snpsum.col, snpsum.row, genofile, clinical, pcs, epidemiological_model, mixedPlot, phenodata, file= "project3_save.RData")

# PART 4 =======================================================================

rm(list = ls())
load("project3_save.RData")
# =============================
# CLEAN GWAS ANNOTATION PIPELINE
# =============================

# Libraries
library(dplyr)
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(gt)

# =============================
# 1. LOAD DATA
# =============================

GWASout <- read.table("GWAA_final.txt", header = FALSE,
                      colClasses = c("character", rep("numeric", 4)))

colnames(GWASout) <- c("SNP", "effect", "stdErr", "tValue", "pValue")

GWASout <- merge(
  GWASout,
  genoMap[, c("SNP", "chr", "position", "A1", "A2")],
  by = "SNP"
)

GWASout$MAF <- snpsum.col[GWASout$SNP, "MAF"]

# =============================
# 2. SELECT TOP HITS
# =============================

alpha <- 0.05
bonf <- alpha / nrow(GWASout)

topHits <- GWASout %>%
  dplyr::filter(pValue < bonf)

# =============================
# 3. BUILD GRANGES
# =============================

library(GenomicRanges)

snps_gr <- GRanges(
  seqnames = paste0("chr", topHits$chr),
  ranges = IRanges(start = topHits$position,
                   end   = topHits$position),
  SNP = topHits$SNP
)

# =============================
# 4. GENE + FUNCTIONAL ANNOTATION
# =============================

library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Genes

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genesGR <- genes(txdb)

# Variant functional annotation
# (UTR, intron, exon, splice, etc)

loc_results <- locateVariants(snps_gr, txdb, AllVariants())

loc_df <- data.frame(
  SNP = mcols(snps_gr)$SNP[as.integer(mcols(loc_results)$QUERYID)],
  LOCATION = mcols(loc_results)$LOCATION,
  GENEID = mcols(loc_results)$GENEID
)

# =============================
# 5. CLASSIFY TYPE OF VARIANT EFFECT
# =============================

# Convert LOCATION into biological impact categories
loc_df$MUTATION_TYPE <- dplyr::case_when(
  loc_df$LOCATION == "coding" ~ "nonsynonymous or coding variant",
  loc_df$LOCATION == "intron" ~ "intronic variant",
  loc_df$LOCATION == "fiveUTR" ~ "5' UTR variant",
  loc_df$LOCATION == "threeUTR" ~ "3' UTR variant",
  loc_df$LOCATION == "spliceSite" ~ "splice site variant",
  loc_df$LOCATION == "promoter" ~ "regulatory promoter variant",
  TRUE ~ "intergenic or unknown"
)

# =============================
# 6. DISTANCE TO NEAREST GENE
# =============================

hits <- distanceToNearest(snps_gr, genesGR)

annot <- data.frame(
  SNP = mcols(snps_gr)$SNP[queryHits(hits)],
  gene_id = genesGR$gene_id[subjectHits(hits)],
  distance = mcols(hits)$distance
)

# =============================
# 7. GENE SYMBOL ANNOTATION
# =============================

symbolMap <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = annot$gene_id,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

annot$GENESYMBOL <- symbolMap$SYMBOL[
  match(annot$gene_id, symbolMap$ENTREZID)
]

# =============================
# 8. MERGE EVERYTHING
# =============================

final_annot <- dplyr::left_join(
  loc_df,
  annot,
  by = "SNP"
)

# Optional merge with GWAS stats
final_annot <- dplyr::left_join(
  final_annot,
  GWASout,
  by = "SNP"
)

# =============================
# 8.5 OUTPUT
# =============================

final_annot <- final_annot %>% arrange(pValue)

write.csv(final_annot, "GWAS_functional_annotation.csv", row.names = FALSE)

# =============================
# 7. VALIDATION CHECK
# =============================

stopifnot(
  all(as.character(seqnames(snps_gr)) == paste0("chr", topHits$chr))
)

# =============================
# 8. DISPLAY TABLE
# =============================

final_annot %>%
  gt() %>%
  tab_header(
    title = "Top GWAS Hits",
    subtitle = "Clean annotation pipeline (no duplicated gene mapping)"
  ) %>%
  fmt_number(columns = c(effect, stdErr, MAF), decimals = 3) %>%
  fmt_scientific(columns = pValue, decimals = 2) %>%
  cols_label(
    SNP = "rsID",
    chr = "Chromosome",
    position = "Position",
    A1 = "Allele 1",
    A2 = "Allele 2",
    MAF = "MAF",
    effect = "Effect Size",
    stdErr = "Std. Error",
    pValue = "P-value",
    distance = "Distance to gene",
    GENESYMBOL = "Gene Symbol"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Part 5 ==========================================================

# BiocManager::install("bigsnpr")
library(bigsnpr)

bedFile <- path.expand("~/Documents/Bio_Statistics/Lab1Files/GWAS_data.bed")

obj <- snp_attach(sub(".bed", ".rds", bedFile))

G   <- obj$genotypes
map <- obj$map

alpha <- 0.05
bonf <- alpha / nrow(GWASout)

sigSNPs <- GWASout %>%
  dplyr::filter(pValue < bonf)
sigSNPs_chr4 <- sigSNPs[sigSNPs$chr == 4 & sigSNPs$chr == 2, ]
sigSNPs <- sigSNPs %>% arrange(chr)
ind <- match(sigSNPs$SNP, map$marker.ID)
ind <- ind[!is.na(ind)]

ld <- snp_cor(
  G,
  ind.col = ind
)

snp_names <- map$marker.ID[ind]

dimnames(ld) <- list(snp_names, snp_names)

ld_mat <- as.matrix(ld)
ld_df <- as.data.frame(as.table(ld_mat))

snp_order <- rownames(ld_mat)

hm_annot <- final_annot %>% arrange(GENESYMBOL)

hm_annot <- hm_annot %>%
  distinct(SNP, .keep_all = TRUE) 

gene_vec <- hm_annot$GENESYMBOL[
  match(snp_order, hm_annot$SNP,)
]

chr_vec <- hm_annot$chr[
  match(snp_order, hm_annot$SNP,)
]

annot_df <- data.frame(
  SNP = snp_order,
  gene = gene_vec
)

ld_df$gene1 <- annot_df$gene[match(ld_df$Var1, annot_df$SNP)]
ld_df$gene2 <- annot_df$gene[match(ld_df$Var2, annot_df$SNP)]

ld_df$block <- ifelse(ld_df$gene1 == ld_df$gene2,
                      ld_df$gene1,
                      "between_genes")
# 
# library(ggplot2)
# 
# ggplot(ld_df, aes(Var1, Var2, fill = Freq)) +
#   geom_tile() +
#   scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90))
# 


library(ComplexHeatmap)
library(circlize)


ha <- HeatmapAnnotation(
  Gene = gene_vec,
  Chromosome = chr_vec,
  col = list(
    Gene = structure(
      c("#70d6ff", "#ff70a6","#ff9770","#ffd670"),
      names = unique(gene_vec)
    ),
    Chromosome = structure(
      c("#fe7f2d", "#619b8a"),
      names = unique(chr_vec)
    )
  ),
  show_annotation_name = TRUE
)


png("heatmap_all.png", width = 1000, height = 1000, res = 120)

hm_all <- Heatmap(
  ld_mat,
  name = "LD (r²)",
  col = colorRamp2(c(0, 0.5, 1), c("#d62828", "#4ea8de", "#a7c957")),
  top_annotation = ha,
  rect_gp = gpar(type = "none"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    # Keep only upper triangle + diagonal
    if(i <= j) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      grid.text(sprintf("%.1f", ld_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    }
  }
)

# 3. Draw the heatmap
draw(hm_all)

# 4. Close the device
dev.off()

library(ggplot2)
library(ggrepel)  # Labels

# 1. MANHATTAN PLOT ================================================


manhattan_data <- GWASout
manhattan_data$CHR <- as.numeric(manhattan_data$chr)
manhattan_data$BP <- manhattan_data$position
manhattan_data$P <- manhattan_data$pValue

chr_info <- manhattan_data %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

manhattan_data <- manhattan_data %>%
  inner_join(chr_info %>% dplyr::select(CHR, offset), by = "CHR") %>%
  mutate(BP_cum = BP + offset)

# Calculate the mean
axis_set <- manhattan_data %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum))

bonferroni <- -log10(0.05 / nrow(GWASout))
suggestive <- -log10(1e-5)


manhattan_data$Neg_logP <- -log10(manhattan_data$P)

top_snps <- manhattan_data %>%
  arrange(P) %>%
  head(5)


# Plot
manhattan_plot <- ggplot(manhattan_data, 
                         aes(x = BP_cum, y = Neg_logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = rep(c("#fbf8cc", "#fde4cf", "#ffcfd2","#f1c0e8","#cfbaf0","#a3c4f3",
                                    "#90dbf4","#8eecf5","#98f5e1","#b9fbc0","#E7F598"), 2)) +
  
  geom_hline(yintercept = bonferroni, linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  geom_hline(yintercept = suggestive, linetype = "dashed", 
             color = "blue", linewidth = 0.6) +
  
  geom_label_repel(data = top_snps,
                   colour = "#f72585",
                   fontface = "bold",
                   aes(label = SNP),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.3,
                   segment.color = "grey50",
                   max.overlaps = 20) +
  
  scale_x_continuous(breaks = axis_set$center, labels = axis_set$CHR) + 

  labs(x = "Chromosome",
       y = expression(-log[10](italic(P))),
       title = "Genome-Wide Association Study",
       subtitle = paste0("Serum Urate levels (n=", nrow(phenodata), ")"),
       caption = paste0("Red line: Bonferroni correction (P=", 
                        signif(0.05/nrow(GWASout), 2), ")")) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 0)
  )

print(manhattan_plot)

ggsave("manhattan_plot.png", manhattan_plot, 
       width = 13, height = 6, dpi = 300)


# 2. QQ PLOT =======================================================

library(ggplot2)

# =========================
# DATA QQ PLOT 
# =========================

p <- GWASout$pValue
p <- p[!is.na(p)]
n <- length(p)

observed <- -log10(sort(p))

expected <- -log10(ppoints(n))

qq_data <- data.frame(
  observed = observed,
  expected = expected
)

# =========================
# LAMBDA GC
# =========================

chisq <- qchisq(1 - p, df = 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

alpha <- seq_len(n)
beta <- n - alpha + 1

ci <- 0.95

upper <- -log10(qbeta((1 + ci)/2, alpha, beta))
lower <- -log10(qbeta((1 - ci)/2, alpha, beta))

qq_data$upper <- upper
qq_data$lower <- lower

qq_data$deviation <- qq_data$observed - qq_data$expected

# =========================
# QQ PLOT
# =========================

qqplot <- ggplot(qq_data, aes(x = expected, y = observed)) +
  
  geom_point(aes(color = deviation),alpha = 0.6, size = 1.8) +
  
  scale_colour_gradient2(mid = "#4361ee", low = "grey70" ,high = "#ff4d6d") +
  
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
            fill = "grey70", alpha = 0.3)+
  
  
  labs(
    x = expression(Expected~~-log[10](italic(P))),
    y = expression(Observed~~-log[10](italic(P))),
    title = "Q-Q Plot",
    subtitle = paste0("lambda GC = ", round(lambda, 3))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank()
  )


ggsave("qq_plot.png", qqplot, 
       width = 13, height = 6, dpi = 300)



# 3. REGIONAL ASSOCIATION PLOT (LOCUS ZOOM) ================================

# Enfocarse en el top hit
GWASout$Neg_logP <- -log10(GWASout$pValue)
top_variant <- GWASout[which.min(GWASout$pValue), ]
locus_chr <- top_variant$chr
locus_pos <- top_variant$position

# Definir ventana de 500kb alrededor del top hit
window <- 500000
locus_data <- GWASout %>%
  filter(chr == locus_chr,
         position >= locus_pos - window,
         position <= locus_pos + window)

region <- GRanges(
  seqnames = paste0("chr", locus_chr),
  ranges = IRanges(start = locus_pos - window,
                   end = locus_pos + window)
)

genes_in_region <- subsetByOverlaps(genesGR, region)

genes_in_region <- as.data.frame(genes_in_region)

genes_in_region$start <- as.numeric(genes_in_region$start)
genes_in_region$end <- as.numeric(genes_in_region$end)


symbolMap <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes_in_region$gene_id,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

genes_in_region$GENESYMBOL <- symbolMap$SYMBOL[
  match(genes_in_region$gene_id, symbolMap$ENTREZID)
]

regional_plot <- ggplot(locus_data, aes(x = position/1e6, y = Neg_logP)) +
  
  # Puntos coloreados por significancia
  geom_point(aes(color = Neg_logP > bonferroni), 
             size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = c("#073b4c", "#ef476f"),
                     labels = c("Not significant", "Genome-wide significant")) +
  
  # Destacar el top SNP
  geom_point(data = top_variant,
             aes(x = position/1e6, y = Neg_logP),
             color = "#06d6a0", size = 4, shape = 18) +
  
  geom_label_repel(data = top_variant,
                   aes(x = position/1e6, y = Neg_logP, label = SNP),
                   color = "#06d6a0", fontface = "bold") +
  
  # Línea de significancia
  geom_hline(yintercept = bonferroni, 
             linetype = "dashed", color = "red") +
  
  # Marcar genes
  geom_segment(data = genes_in_region,
               aes(x = start/1e6, xend = end/1e6,
                   y = -0.5, yend = -0.5),
               color = "#c1121f", linewidth = 3) +
  
  geom_text(data = genes_in_region,
            aes(x = (start + end)/(2*1e6), y = -1, label = GENESYMBOL),
            size = 3, angle = 0, color = "#c1121f", fontface = "italic") +
  
  labs(x = paste0("Position on chromosome ", locus_chr, " (Mb)"),
       y = expression(-log[10](italic(P))),
       title = paste0("Regional Association Plot - chr", locus_chr),
       subtitle = paste0("±500kb around ", top_variant$SNP),
       color = "") +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(regional_plot)

ggsave("regional_plot.png", regional_plot, 
       width = 20, height = 7, dpi = 300)


# 4. EFFECT SIZE PLOT =======================================================

top20 <- GWASout %>%
  filter(pValue < 0.05/nrow(GWASout)) %>%
  arrange(pValue) %>%
  head(20)

if(nrow(top20) > 0) {
  
  top20$SNP_label <- paste0(top20$SNP, " (", top20$A2, ")")
  top20$SNP_label <- factor(top20$SNP_label, 
                            levels = top20$SNP_label[order(top20$effect)])
  
  effect_plot <- ggplot(top20, aes(x = SNP_label, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    
    geom_errorbar(aes(ymin = effect - qnorm(0.975)*stdErr,
                      ymax = effect + qnorm(0.975)*stdErr),
                  width = 0.3, color = "grey30", linewidth = 0.8) +
    
    geom_point(aes(color = effect > 0), size = 4) +
    
    scale_color_manual(values = c("#D62828", "#0077B6"),
                       labels = c("Risk", "Protective")) +
    
    coord_flip() +
    
    labs(x = "Variant (Effect Allele)",
         y = "Effect Size (β) ± 95% CI",
         title = "Effect Sizes of Genome-Wide Significant Variants",
         color = "Direction") +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  print(effect_plot)
  
  ggsave("effect_size_plot.png", effect_plot, 
         width = 10, height = 8, dpi = 300)
}



# 5. SUMMARY TABLE (PRETTY) ===================================================

summary_table_clean <- final_annot %>%
  distinct(SNP, .keep_all = TRUE) 

summary_table <- summary_table_clean %>%
  mutate(
    `95% CI` = paste0("[", 
                      round(effect - 1.96*stdErr, 3), ", ",
                      round(effect + 1.96*stdErr, 3), "]")
  ) %>%
  dplyr::select(SNP, chr, position, GENESYMBOL, 
         A1, A2, MAF, effect, `95% CI`, pValue, LOCATION, MUTATION_TYPE) %>%
  gt() %>%
  tab_header(
    title = md("**Top GWAS Hits for Serum Urate**"),
    subtitle = md("*Genome-wide significant associations*")
  ) %>%
  fmt_number(columns = c(MAF, effect), decimals = 3) %>%
  fmt_scientific(columns = pValue, decimals = 2) %>%
  fmt_number(columns = position, decimals = 0, sep_mark = ",") %>%
  cols_label(
    SNP = md("**rsID**"),
    chr = md("**Chr**"),
    position = md("**Position**"),
    GENESYMBOL = md("**Gene**"),
    LOCATION = md("**Location**"),
    MUTATION_TYPE = md("**Mutation Type**"),
    A1 = md("**A1**"),
    A2 = md("**A2**"),
    MAF = md("**MAF**"),
    effect = md("**β**"),
    `95% CI` = md("**95% CI**"),
    pValue = md("***P*-value**")
  ) %>%
  data_color(
    columns = pValue,
    colors = scales::col_numeric(
      palette = c("#560bad", "#b5179e", "#f72585"),
      domain = NULL,
      reverse = TRUE
    )
    
  ) %>%
  data_color(
    columns = LOCATION,
    colors = scales::col_factor(
      palette = c(
        "intron" = "#560bad",
        "intergenic" = "#b5179e",
        "threeUTR" = "#f72585",
        "fiveUTR" = "#7209b7",
        "promoter" = "#f15bb5",
        "coding" = "#00bbf9",
        "spliceSite" = "#00f5d4"
      ),
      domain = NULL
    )
  ) %>%
  data_color(
    columns = GENESYMBOL,
    colors = scales::col_factor(
      palette = c("EIF2B4" = "#70d6ff", "GCKR" = "#ff70a6","GTF3C2" = "#ff9770","SLC2A9" = "#ffd670"),
      domain = NULL
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = SNP)
  ) %>%
  tab_options(
    table.font.size = 11,
    heading.align = "center",
    column_labels.font.weight = "bold"
  ) %>%
  tab_source_note(
    source_note = paste0("Bonferroni threshold: P < ", 
                         signif(0.05/nrow(GWASout), 2))
  )

print(summary_table)

gtsave(summary_table, "gwas_summary_table.html")
gtsave(summary_table, "gwas_summary_table.png")



# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

clinical$sex <- factor(
  clinical$sex,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

histp <- clinical %>%
  ggplot(aes(x = Serum_Urate_mgdL, fill = sex)) +
  
  theme(
    legend.position = "top"
  )+
  
  geom_histogram(
    color = "#e9ecef",
    alpha = 0.6,
    position = "identity",
    bins = 30
  ) +
  
  geom_vline(
    xintercept = 7,
    color = "#4895ef",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  
  geom_vline(
    xintercept = 6,
    color = "#ff0054",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  
  scale_fill_manual(values = c("#780000", "#390099")) +
  theme_ipsum() +
  labs(fill = "")+

  annotate("text", x = 7, y = 50,
         label = "Men", color = "#001524", angle = 90)+

  annotate("text", x = 6, y = 50,
         label = "Women", color = "#001524", angle = 90, fontface = "bold")

print(histp)


ggsave("histogram.png", histp, 
       width = 10, height = 8, dpi = 300)


##==============================================================================
library(ggplot2)
library(dplyr)

# =========================================================
# MODEL
# =========================================================

model_summary <- summary(epidemiological_model)

coef_table <- as.data.frame(model_summary$coefficients)
coef_table$variable <- rownames(coef_table)

# =========================================================
# FILTER
# =========================================================

sig_coef <- coef_table %>%
  filter(
    `Pr(>|t|)` < 0.05,
    variable != "(Intercept)"
  )

# =========================================================
# RENAME
# =========================================================

sig_coef$variable <- c(
  "Be female",
  "Age",
  "Alcohol",
  "Use Diuretics",
  "Exercise",
  "BMI*Exercise"
)

# =========================================================
# SCALE
# =========================================================

min_abs <- min(abs(sig_coef$Estimate))
max_abs <- max(abs(sig_coef$Estimate))

sig_coef <- sig_coef %>%
  mutate(
    
    # Escalado visual
    Estimate_scaled =
      sign(Estimate) *
      (
        0.15 +
          0.85 *
          (
            (abs(Estimate) - min_abs) /
              (max_abs - min_abs)
          )
      ),
    
    # Dirección para colores
    direction = ifelse(
      Estimate > 0,
      "Positive",
      "Negative"
    )
  )

# =========================================================
# SORT
# =========================================================

sig_coef$variable <- factor(
  sig_coef$variable,
  levels = sig_coef$variable[
    order(sig_coef$Estimate_scaled)
  ]
)

sig_coef <- sig_coef %>%
  mutate(
    
    significance = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01  ~ "**",
      `Pr(>|t|)` < 0.05  ~ "*",
      TRUE ~ ""
    )
    
  )
# =========================================================
# GRAPH
# =========================================================

estimates_plot <- ggplot(
  sig_coef,
  aes(
    x = variable,
    y = Estimate_scaled,
    fill = direction
  )
) +
  
  # Barras
  geom_col(
    width = 0.7
  ) +
  
  # Línea central
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # P-values
  geom_text(
    aes(
      label = paste0(
        significance,
        "\n",
        "p=",
        signif(`Pr(>|t|)`, 2)
      )
    ),
    
    nudge_y = ifelse(
      sig_coef$Estimate > 0,
      0.08,
      -0.08
    ),
    
    size = 4.5,
    fontface = "bold"
  )+
  
  # Voltear ejes
  coord_flip(
    clip = "off"
  ) +
  
  # Colores
  scale_fill_manual(
    values = c(
      "Positive" = "#3b82f6",
      "Negative" = "#ef4444"
    )
  ) +
  
  # Límites visuales
  expand_limits(
    y = c(-1.2, 1.2)
  ) +
  
  # Tema
  theme_minimal(
    base_size = 18
  ) +
  
  theme(
    
    legend.position = "none",
    
    panel.grid.major.y =
      element_blank(),
    
    axis.title.y =
      element_blank(),
    
    plot.margin =
      margin(20, 80, 20, 20),
    
    plot.title =
      element_text(
        face = "bold"
      )
  ) +
  
  # Etiquetas
  labs(
    title = "Covariates",
    y = "(Min-Max)Scaled Effect Size"
  )

save(bonferroni, hm_all, manhattan_plot, qqplot, regional_plot, effect_plot, summary_table, estimates_plot, file= "project4_save.RData")

# ==============================================================================
