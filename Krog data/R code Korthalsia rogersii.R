# The script used for this analysis was adapted from existing resources, including tutorials such as:(a) https://bookdown.org/hhwagner1/LandGenCourse_book/, (b)https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/ 
library(adegenet) # for data converstion, manipulation, Descriptive Statistics and DAPC
library(poppr) # Descriptive Statistics
library(dplyr)
library(hierfstat) # Descriptive Statistics
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
# read genalex data as table (.txt)
genalex_data <- read.table("/K.rogersii 1.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
# creat a funtion for grouping locus 
combine_alleles <- function(locus1, locus2) {
  return(paste(locus1, locus2, sep = "/"))}
genalex_data$kle16 <- mapply(combine_alleles, genalex_data$Kle16_1, genalex_data$Kle16_2)
genalex_data$kle18 <- mapply(combine_alleles, genalex_data$Kle18_1, genalex_data$Kle18_2)
genalex_data$kle23 <- mapply(combine_alleles, genalex_data$Kle23_1, genalex_data$Kle23_2)
genalex_data$kle30 <- mapply(combine_alleles, genalex_data$Kle30_1, genalex_data$Kle30_2)
genalex_data$kle3 <- mapply(combine_alleles, genalex_data$Kle3_1, genalex_data$Kle3_2)
genalex_data$kle13 <- mapply(combine_alleles, genalex_data$Kle13_1, genalex_data$Kle13_2)
genalex_data$klest24 <- mapply(combine_alleles, genalex_data$klest24, genalex_data$klest24)
b<- data.frame(genalex_data)
View(b)
a <- subset(b, select = c(Sample, Pop,Kle16,Kle18,Kle23,Kle30,Kle3,Kle13,Klest24))
krog_gen <- df2genind(a[, c("Kle16","Kle18","Kle23","Kle30","Kle3","Kle13","Klest24")],
                        sep = "/", NA.char = "0",pop = a$Pop,ind.names = a$Sample, ploidy = 2,ncode = 3L)  
saveRDS(krog_gen, file = "krog_gen.rds")
krog_gen <- readRDS("krog_gen.rds")

# Basic popgen analysis 
# Data Filtering 
# Missing data per loci 
locmiss_krog = propTyped(krog_gen, by = "loc")
locmiss_krog[which(locmiss_krog < 0.80)]
krog_gen = missingno(krog_gen, type = "loci", cutoff = 0.20) #Remove microsatellite loci with > 20% missing data
# Missing data per individuals
indmiss_krog= propTyped(krog_gen, by = "ind")
indmiss_krog[ which(indmiss_krog < 0.80) ] # print individuals with < 80% complete genotypes
krog_gen = missingno(krog_gen, type = "geno", cutoff = 0.20)
# Check unique genotypes 
mlg(krog_gen)
dups_krog = mlg.id(krog_gen)
for (i in dups_krog){ # for each element in the list object
  if (length(dups_krog[i]) > 1){ # if the length is greater than 1
    print(i) # print individuals that are duplicates
  }
  }
krog_dups <- c("Gkr3","Bkr8","Gkr22","Ckr12","Dkr10") # duplicate genotype identified
krog_Nodups = indNames(krog_gen)[! indNames(krog_gen) %in% krog_dups]
krog_gen = krog_gen[krog_Nodups, ]
mlg(krog_gen)
isPoly(krog_gen) %>% summary # check polymorphism

# Check for deviations from Hardy-Weinberg equilibrium (HWE)
# HWE globally across individuals
round(pegas::hw.test(krog_gen, B = 1000), digits = 3)
# HWE of each locus in each population
# Chi-squared test: p-value
HWE.test <- data.frame(sapply(seppop(krog_gen), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
  round(HWE.test.chisq,3)}

# Monte Carlo: p-value
HWE.test <- data.frame(sapply(seppop(krog_gen), 
                              function(ls) pegas::hw.test(ls, B=1000)[,4]))
HWE.test.MC <- t(data.matrix(HWE.test))
{cat("MC permuation test (p-values):", "\n")
  round(HWE.test.MC,3)}
##proportion of populations out of Hardy-Weinberg equilibrium (HWE), conservative cut-off of alpha = 0.05 for each test

alpha=0.05
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                   MC=apply(HWE.test.MC<alpha, 2, mean))
Prop.loci.out.of.HWE

#proportion of populations out of Hardy-Weinberg equilibrium (HWE)
Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean))
Prop.pops.out.of.HWE             


#Summary statistics
krog_gen
table(krog_gen$loc.fac) # number of alleles per locus 
summary(krog_gen$pop) # pop size 
private_alleles(krog_gen) %>% apply(MARGIN = 1, FUN = sum) #number of private alleles 
allelic.richness(genind2hierfstat(krog_gen))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3) #allelic richness per population
basic_krog = basic.stats(krog_gen, diploid = TRUE) #Calculate basic stats using hierfstat
# Oberved heteozygosity per population 
Ho_krog = apply(basic_krog$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Ho_krog
# Expected heteozygosity per population 
He_krog = apply(basic_krog$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
He_krog
# Inbreeding coefficient (FIS)
apply(basic_krog$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
# pairwise FST 
krog_fst = genet.dist(krog_gen, method = "WC84") %>% round(digits = 3)
krog_fst 


#discriminant analysis of principal components (DAPC) 
# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
x = tab(krog_gen, NA.method = "mean")
crossval = xvalDapc(x, krog_gen$pop, result = "groupMean", xval.plot = TRUE)
# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
crossval$`Number of PCs Achieving Highest Mean Success`
crossval$`Number of PCs Achieving Lowest MSE`
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)
# Run a DAPC using site IDs as priors
dapc1 = dapc(krog_gen, krog_gen$pop, n.pca = numPCs, n.da = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))
# Create a data frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(krog_gen)

# Add a column with the site IDs
ind_coords$Site = krog_gen$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(krog_gen), "Set2")
# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) )

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("Korthalsia rogersii DAPC")+
  # custom theme
  ggtheme
