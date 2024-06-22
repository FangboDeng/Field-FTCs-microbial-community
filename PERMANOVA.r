##bacterial communites analysis as an example

library(vegan)
###PERMANOVA 
map<-read.table("otu_mapping.txt",header=F,sep="\t")#mapping文件里的name加#
FTC<-as.vector(map[,2])
stover<-as.vector(map[,3])
fertilizer<-as.vector(map[,4])

#####Three factors  
set.seed(666)
BC <- read.table(file="Bray Curtis distance-matrix.tsv")
BCdist <- as.dist(BC)
adonis1=adonis2(BCdist~FTC*stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

set.seed(666)
unifrac <- read.table(file="Weighted Unifrac distance-matrix.tsv")
unifracdist <- as.dist(unifrac)
adonis2=adonis2(unifracdist~FTC*stover*fertilizer,permutations = 9999) #用adonis互作
adonis2

####Two factors
##Autumn
map<-read.table("otu_mapping BF.txt",header=F,sep="\t")#mapping文件里的name加#
stover<-as.vector(map[,3])
fertilizer<-as.vector(map[,4])

BC <- read.table(file="BF Bray Curtis distance-matrix.tsv")
BCdist <- as.dist(BC)
adonis1=adonis2(BCdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

unifrac <- read.table(file="BF Weithted Unifrac distance-matrix.tsv")
unifracdist <- as.dist(unifrac)
adonis1=adonis2(unifracdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

##Winter
map<-read.table("otu_mapping DF.txt",header=F,sep="\t")#mapping文件里的name加#
stover<-as.vector(map[,3])
fertilizer<-as.vector(map[,4])

BC <- read.table(file="DF Bray Curtis distance-matrix.tsv")
BCdist <- as.dist(BC)
adonis1=adonis2(BCdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

unifrac <- read.table(file="DF Weithed Unifrac distance-matrix.tsv")
unifracdist <- as.dist(unifrac)
adonis1=adonis2(unifracdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

##Spring
map<-read.table("otu_mapping TF.txt",header=F,sep="\t")#mapping文件里的name加#
stover<-as.vector(map[,3])
fertilizer<-as.vector(map[,4])

BC <- read.table(file="TF Bray Curtis distance-matrix.tsv")
BCdist <- as.dist(BC)
adonis1=adonis2(BCdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1

unifrac <- read.table(file="TF Weithed Unifrac distance-matrix.tsv")
unifracdist <- as.dist(unifrac)
adonis1=adonis2(unifracdist~stover*fertilizer,permutations = 9999) #用adonis互作
adonis1
