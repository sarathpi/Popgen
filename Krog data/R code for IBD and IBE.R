#IBD IBE
#requied packages 
library(vegan) # for mantel and partial mantel test 
library(geosphere)# for GPS to geographical distance matrix 
library(geodist) # for GPS to geographical distance matrix 
library(gstudio) # for importing genetic data (genepop) format ape also can be used for this 
require(devtools)
library(ape) #for calculating pairwise genetic distance
library (raster)
library (terra) #for importing raster layer 
a <- read_population("../krog_genepop.txt", type = "genepop")
Genetic.dist <- dist.gene(a,method = "pairwise",pairwise.deletion = FALSE,variance = FALSE)
Genetic.dist
class(Genetic.dist)
Genetic.dist.matrix <- Genetic.dist
#import geographical coordinates for all individuals 
gdis <-read.csv(".../GPS.csv",row.names = 1)
#calculate geographical distance matrix
gdis.dist.matrix <- as.dist(geodist(gdis,measure = "haversine"))
#import and extract raster values for all environmental layers 
raster <- rast(c(
  "../bio_19.tif",
  "../bio_1.tif",
  "../bio_2.tif",#import all 19 bioclimatic variables))

#Converitng gps coords to spatial points
points <- vect(gdis, geom = c("lon", "lat"), crs = raster)

#extract env values for each coord
values <- extract(raster, points)
values

#select the data subset
env2 <- values[2:20]
#scale and center the matrix
envss <- scale(env2, center = TRUE)

#env distance matrix
envdist <- dist(envss,method = "euclidean")
#to matrix format 
env.dist.matrix <- as.dist(as.matrix(envdist))
#names 
names(Genetic.dist.matrix) <-names(gdis.dist.matrix)<-names(env.dist.matrix)

#mantel test
g.gd.matel <- mantel(gdis.dist.matrix,Genetic.dist.matrix,method = "pearson", permutations = 1000, strata = NULL, na.rm = TRUE, parallel = getOption("mc.cores"))

g.env.matel <- mantel(Genetic.dist.matrix,env.dist.matrix,method = "pearson", permutations = 1000, strata = NULL, na.rm = TRUE, parallel = getOption("mc.cores"))

#partial mantel test

g.gd.env.partialmatel <- mantel.partial(Genetic.dist.matrix,gdis.dist.matrix, env.dist.matrix, method = "pearson", permutations = 100000,strata = NULL, na.rm = TRUE, parallel = getOption("mc.cores"))

g.env.gd.partialmatel <- mantel.partial(Genetic.dist.matrix, env.dist.matrix,gdis.dist.matrix, method = "pearson", permutations = 100000,strata = NULL, na.rm = TRUE, parallel = getOption("mc.cores"))
