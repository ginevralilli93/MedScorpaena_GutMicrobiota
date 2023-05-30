# Diversity analyses (Alpha and Beta) as reported in the article Lilli et al., .

# The phyloseq object called in this script has been filtred previously for taxa_sums=0, mithocondria, cloroplasts and Archea, samples_sums<10000
# Comment: The taxaused are classified with the RDP classifier using the silva v138 database
# Species dataset (only BO) used for testing Species effect. Analyzed at the ASV level. 
# Extended dataset used for all the Geo analyses for each species separately --> datset including data from 3 studies = this study, Kormas et al., 2022, Ruiz-Rodriguez et al., 2020
# Extended dataset is analyzed only at the Genus level. 

# Author: Ginevra Lilli 


# Libraries--------------------------------------------------------------------------
#Loading the packages that we will use
pkgs <- c("tidyverse", "phyloseq", "ggpubr", "ggplot2", 
          "vegan", "reshape2", "DESeq2","microbiome","colorspace",
          "picante","RColorBrewer","ggpubr","zCompositions","ggord",
          "forecast","rstatix","multcomp",
          "MicEco")
lapply(pkgs, require, character.only = TRUE)


# Generate the palettes--------------------------------------------------------------
species_pal2<-c("royalblue4","orange2","yellow2")
mpa.palette.sn<-c("#FFCC00","grey60","grey2")
mpa.palette.sp.ss<-c("grey2","tan1")
theme_set(theme_minimal())
zone.palette<-sample(colors(),10)
# Functions------------------------------------------------------------------------

# Phy object rarefaction #
tran.rare<-function(phy.obj, sample.size){
  rare<-rarefy_even_depth(phy.obj,rngseed = 123,sample.size=sample.size, replace=FALSE)
  #select only the asv that have colSums !=0. In this way the next steps can be performed
  colsum.nozero<-(colSums(otu_table(rare), na.rm=T) !=0)
  nozero<-otu_table(phy.obj)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  filt.obj<-phyloseq(otu_table(nozero, taxa_are_rows = FALSE),
                     tax_table(tax_table(rare)), taxa_names(taxa_names(rare)),
                     sample_data(sample_data(rare)),phy_tree(phy_tree(rare)))
}


#create the alpha div table 
alpha.div<-function(phy,metadata){ #phy: phyloseq object rarefied
  alpha_div<- phy %>% 
    estimate_richness(split= T, measures = c("Shannon","Observed")) %>%
    mutate(X.SampleID = sample_names(phy)) # add sample IDs
  colnames(alpha_div)[3]<-"FishID"
  #estimate richness by phylogeny: Faith's index - calculate alpha diversity 
  tree<-phy_tree(phy)
  asv.table<-otu_table(phy)
  faith.index<-pd(asv.table,tree,include.root = FALSE)
  faith.index$FishID<-rownames(faith.index)
  
  #add the next column with the alpha diversity indexes in the metadata table. 
  alpha_div <- cbind(metadata, # add metadata
                     alpha_div,
                     faith.index,
                     by = "FishID")
  alpha_div<-alpha_div[,-c(63:64)] #remove "by" and "FishID" since there are several columns called like that
}

# Generate a plot for the alpha diversity values
alpha.plot<-function(alpha_div,palette,method, points, color.points){
  alpha_div %>% 
    ggplot(aes(as.character(MPA),Shannon))+
    geom_boxplot(data=alpha_div,  aes(alpha=0.5), fill=palette) + 
    geom_point(data= points, aes(as.character(MPA),Shannon), color ={{color.points}}, alpha =0.3, size=4)+
    stat_compare_means(method = {{method}})+ 
    #scale_x_discrete(labels=c("1" = "Banyuls-sur-mer", "2" = "Carry-Le-Rouet", "3"="Bonifacio", "4"="Gyaros island"))+
    #facet_wrap(~Species)+
    ylab("Shannon's index")+ xlab("")+ylim(0,8)+
    # stat_summary(fun.y="mean", color ="red", shape= 20)+
    theme(legend.key.size = unit(4, "mm"), # resize legend so it fits in the figure
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 7),
          axis.text.y = element_text(size=15),
          axis.title.y = element_text(size=20,face="bold"),
          axis.text.x = element_text(angle = 45,hjust=0.9, size = 18, face = "italic"))
  
}

#Centered log-ratio trasformation
tran.clr<-function(phy.obj){
  
  #agglomerate at the genus level
  phy.obj.glom<-tax_glom(phy.obj,taxrank = "Genus",NArm = FALSE)
  #select only the asv that have colSums !=0. In this way the next steps can be performed
  colsum.nozero<-(colSums(otu_table(phy.obj.glom), na.rm=T) !=0)
  nozero<-otu_table(phy.obj.glom)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #transform the 0s in probability vectors
  asv.prob<-zCompositions::cmultRepl(nozero, method = 'CZM', delta = 0.5, output = 'p-counts')
  
  #create an alternative ps object with the asv table with probabilities 
  prob<-phyloseq(otu_table(asv.prob, taxa_are_rows = FALSE),
                 tax_table(tax_table(phy.obj.glom)), taxa_names(taxa_names(phy.obj.glom)),
                 sample_data(sample_data(phy.obj.glom)),phy_tree(phy_tree(phy.obj.glom)))
  
  clr <- microbiome::transform(prob, "clr")  
  
}

#Centered log-ratio trasformation - agglomerated to genus 
tran.clr.not.glom<-function(phy.obj){
  
  #agglomerate at the genus level
  #phy.obj.glom<-tax_glom(phy.obj,taxrank = "Genus")
  #select only the asv that have colSums !=0. In this way the next steps can be performed
  colsum.nozero<-(colSums(otu_table(phy.obj), na.rm=T) !=0)
  nozero<-otu_table(phy.obj)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #transform the 0s in probability vectors
  asv.prob<-zCompositions::cmultRepl(nozero, method = 'CZM', delta = 0.5, output = 'p-counts')
  
  #create an alternative ps object with the asv table with probabilities 
  prob<-phyloseq(otu_table(asv.prob, taxa_are_rows = FALSE),
                 tax_table(tax_table(phy.obj)), taxa_names(taxa_names(phy.obj)),
                 sample_data(sample_data(phy.obj)),phy_tree(phy_tree(phy.obj)))

  clr <- microbiome::transform(prob, "clr")  
  
}

# Ordination plots - PCA - Species
pcaSpecies<-function(phy.clr,palette){
  pca.clr <- ordinate(phy.clr, "RDA") #Performs redundancy analysis, or optionally principal components analysis, via rda
  sample_data(phy.clr)$Species<-as.character(sample_data(phy.clr)$Species)
  plot_ordination(phy.clr,pca.clr, color="Species", shape = "Zone")+
    geom_point(size=4, alpha=1, aes(color=as.character(Species)))+
    scale_color_manual(values = palette)+
    #  facet_wrap(~Species)+
    #stat_ellipse(geom = "polygon",  alpha=0.09,size=1)+
    theme(axis.text.x = element_text(size=15))+
    theme(axis.text.y = element_text(size=15))+
    theme(axis.title.x  = element_text(size=15))+
    theme(axis.title.y  = element_text(size=15))+
    theme(legend.text = element_text(face = "italic", size= 15))+
    theme(legend.title = element_text( size= 15))
}

# Ordination plots - PCA - MPA
pcaMPA<-function(phy.clr,palette){
  pca.clr <- ordinate(phy.clr, "RDA") #Performs redundancy analysis, or optionally principal components analysis, via rda
  sample_data(phy.clr)$MPA<-as.character(sample_data(phy.clr)$MPA)
  plot_ordination(phy.clr,pca.clr, color="MPA")+
    geom_point(size=4, alpha=1, aes(color=as.character(MPA) ))+
    scale_color_manual(values = palette)+
    #  facet_wrap(~Species)+
    stat_ellipse(geom = "polygon",  alpha=0.09,size=1)+
    theme(axis.text.x = element_text(size=15))+
    theme(axis.text.y = element_text(size=15))+
    theme(axis.title.x  = element_text(size=15))+
    theme(axis.title.y  = element_text(size=15))+
    theme(legend.text = element_text(face = "italic", size= 15))+
    theme(legend.title = element_text( size= 15))
}

# Ordination plots - PCA - MPA
pcaZone<-function(phy.clr,palette){
  pca.clr <- ordinate(phy.clr, "RDA") #Performs redundancy analysis, or optionally principal components analysis, via rda
  sample_data(phy.clr)$Zone<-as.character(sample_data(phy.clr)$Zone)
  plot_ordination(phy.clr,pca.clr, color="Zone")+
    geom_point(size=4, alpha=1, aes(color=as.character(Zone) ))+
    scale_color_manual(values = palette)+
    #  facet_wrap(~Species)+
    #stat_ellipse(geom = "polygon",  alpha=0.09,size=1)+
    theme(axis.text.x = element_text(size=15))+
    theme(axis.text.y = element_text(size=15))+
    theme(axis.title.x  = element_text(size=15))+
    theme(axis.title.y  = element_text(size=15))+
    theme(legend.text = element_text(face = "italic", size= 15))+
    theme(legend.title = element_text( size= 15))
}




# Welch MANOVA
Welch.Manova<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  WDS<-WdS.test(dm = matrix, f = df[,var], nrep =999) 
}

#PERMANOVA
permanova<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  permanova<-adonis(matrix~df[,var])
}


# PERMDISP - betadisperser

permdisp<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df$MPA<-as.factor(df$MPA) 
  dispr<-betadisper(matrix,df[,var]) #betadisperser test 
  permutest(dispr) 
}

# PAIRWISE ADONIS 
pairwise.adonis<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  pairwise<-pairwiseAdonis::pairwise.adonis(matrix, factors=df[,var]) 
}

# Create a function to remove the taxa with 0 counts in the new phy obj

nozero_taxa<-function(phy_rare){
  colsum.nozero<-(colSums(otu_table(phy_rare), na.rm=T) !=0)
  nozero1<-otu_table(phy_rare)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  phy_rare2<-phyloseq(otu_table(nozero1, taxa_are_rows = FALSE),
                      tax_table(tax_table(phy_rare)), taxa_names(taxa_names(phy_rare)),
                      sample_data(sample_data(phy_rare)),phy_tree(phy_tree(phy_rare)))
  
}

# Gather the data

# Obtain the dada2+phyloseq object with additional samples from the 2 studies

# This is the Extended datset

sf<- get(load("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/merged_object/Rdata_obj/sf.filtered.object.with.tree.RData"))
meta.zone.location<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Data/metadata/125_scorpionfish_zone_by_loc.csv") 
meta.zone.location<-meta.zone.location[,-1]
rownames(meta.zone.location)<-meta.zone.location$FishID
sf<-phyloseq(otu_table(otu_table(sf), taxa_are_rows = FALSE),
             tax_table(tax_table(sf)), taxa_names(taxa_names(sf)),
             sample_data(meta.zone.location),phy_tree(phy_tree(sf)))


# 1 - Results chapter: "The effect of host's phylogeny on the scorpionfish gut microbiota" ####

# Analysis performed at the ASV level

# Subset the phyloseq object for the locations
species.dataset.BO<-subset_samples(sf, MPA %in% 3) # 3 is the identifier of the BO region
species.dataset.GY<-subset_samples(sf, MPA %in% 4) # 4 is the identifier of the GY region

# Rarefy 
species.rare.BO<-tran.rare(species.dataset.BO, sample.size=20000)
species.rare.GY<-tran.rare(species.dataset.GY, sample.size=18025) # Keep 18025 because one of the external samples has a lower number of reads than 20000

#Obtain the metadata from each subset 
species.meta.BO<-sample_data(species.rare.BO) %>% data.frame()
species.meta.GY<-sample_data(species.rare.GY) %>% data.frame()

## 1.1. Size difference----------------------------------------------------------------------------

#check the range of sizes
ggplot(species.meta.BO,aes(Species,Length, fill=Species))+
  geom_boxplot()+
  geom_point()+
  stat_compare_means()+
  scale_fill_manual(values = species_pal2)+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        legend.text = element_text(size=25, face="italic"))


# Test the difference of sizes in the three species

shapiro.test(species.meta.BO$Length) #  normally distributed

res.aov<-aov(Length~Species, data=species.meta.BO)
summary.aov(res.aov) # F value 19.86; P value 7.77e-07
TukeyHSD(res.aov)

# $Species
# diff       lwr       upr     p adj
# Sporcus-Snotata 6.578182  3.582606  9.573757 0.0000100
# Sscrofa-Snotata 8.668182  5.050681 12.285683 0.0000020
# Sscrofa-Sporcus 2.090000 -1.007843  5.187843 0.2410732

## 1.2. Alpha diversity ----------------------------------------------------------------------------

# Generate the alpha diversity table
species.alpha.div.BO<-alpha.div(species.rare.BO,species.meta.BO)
species.alpha.div.GY<-alpha.div(species.rare.GY,species.meta.GY)

# Test the alpha diversity

# BO

# Shannon 
shapiro.test(species.alpha.div.BO$Shannon) # p-value = 0.00734
kruskal.test(Shannon~Species, data = species.alpha.div.BO)## p-value = 0.4568
boxplot(Shannon~Species, data = species.alpha.div.BO)

# Observed 
shapiro.test(species.alpha.div.BO$Observed) # p-value = 3.141e-05
kruskal.test(Observed~Species, data = species.alpha.div.BO)##p-value = 0.7099

# Faith
shapiro.test(species.alpha.div.BO$PD) # p-value = 8.86e-06
kruskal.test(PD~Species, data = species.alpha.div.BO)## p-value = 0.6355



# GY

# Shannon 
shapiro.test(species.alpha.div.GY$Shannon) # p-value = 0.04274
kruskal.test(Shannon~Species, data = species.alpha.div.GY)## p-value = 0.3865

# Observed 
shapiro.test(species.alpha.div.GY$Observed) # p-value = 0.688
summary(aov(Observed~Species, data = species.alpha.div.GY))## p-value = 0.743

# Faith
shapiro.test(species.alpha.div.GY$PD) # p-value = 0.7325
summary(aov(PD~Species, data = species.alpha.div.GY))## p-value = 0.613


# Repeat the tests with the ASVs agglomerate to the Genus level
# Do the results of the tests change? 

species.rare.BO.glom<-tax_glom(species.rare.BO, taxrank = "Genus", NArm = FALSE)
species.rare.GY.glom<-tax_glom(species.rare.GY, taxrank = "Genus", NArm = FALSE)

alpha.BO.glom<-alpha.div(species.rare.BO.glom,species.meta.BO)
alpha.GY.glom<-alpha.div(species.rare.GY.glom,species.meta.GY)

#BO 
# Shannon 
shapiro.test(alpha.BO.glom$Shannon) # p-value = 0.0001
kruskal.test(Shannon~Species, data = alpha.BO.glom)## p-value = 0.3419
boxplot(Shannon~Species, data = alpha.BO.glom)


# BOXPLOT

# highlight point from study Rodriguez- Ruiz in green


my_comparisons<-list(c("Snotata","Sporcus"),c("Snotata","Sscrofa"),c("Sporcus","Sscrofa"))

alpha.BO.glom%>% group_by(Species)%>%summarise(freq=n())

ggp.shannon.SpeciesBO.glom<-
  alpha.BO.glom %>% 
  ggplot(aes(Species,Shannon))+
  geom_boxplot(data=alpha.BO.glom, fill= species_pal2, aes(alpha=0.5))+
  geom_point()+
  stat_compare_means(comparisons = my_comparisons,(aes(label = ..p.signif..)))+
  ylab("Shannon's index")+ xlab("")+ylim(0,8)+
  # stat_summary(fun.y="mean", color ="red", shape= 20)+
  theme(legend.key.size = unit(4, "mm"), # resize legend so it fits in the figure
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45,hjust=0.9, size = 18, face = "italic"))




# Observed 
shapiro.test(alpha.BO.glom$Observed) # p-value = 1.012e-05
kruskal.test(Observed~Species, data = alpha.BO.glom)##p-value = 0.6475

# Faith
shapiro.test(alpha.BO.glom$PD) # p-value = 0.0001
kruskal.test(PD~Species, data = alpha.BO.glom)## p-value = 0.6345

# GY
# Shannon 
shapiro.test(alpha.GY.glom$Shannon) # p-value = 0.7385
summary(aov(Shannon~Species, data = alpha.GY.glom))## p-value = 0.232
boxplot(Shannon~Species, data = alpha.GY.glom)

# Observed 
shapiro.test(alpha.GY.glom$Observed) # p-value =0.6652
summary(aov(Observed~Species, data = alpha.GY.glom))##p-value = 0.678

# Faith
shapiro.test(alpha.GY.glom$PD) # p-value = 0.773
summary(aov(PD~Species, data = alpha.GY.glom))## p-value = 0.66

# 1.2.1 - Testing the Alpha diveristy on the Extended dataset--------------------


# Rarefy 
sf.rare<-tran.rare(sf, sample.size=18025) # Keep 18025 because one of the external samples has a lower number of reads than 20000

#Agglomerate to the genus level the three objects whit the external data 
sf.rare.glom<-tax_glom(sf.rare, taxrank = "Genus", NArm = FALSE)

#Obtain the metadata from each subset 
sf.meta<-sample_data(sf.rare) %>% data.frame()

# TABLE
sf.alpha.div<-alpha.div(sf.rare.glom,sf.meta)

# Shannon 
shapiro.test(sf.alpha.div$Shannon) # p-value = 5.624e-08
kruskal.test(Shannon~Species, data = sf.alpha.div)##  significant p-value = 0.007
FSA::dunnTest(Shannon~Species, data = sf.alpha.div) ## siginf between S.scrofa and S.notata p value = 0.013

# Obsrved 
shapiro.test(sf.alpha.div$Observed) # p-value =9.044e-12
kruskal.test(Observed~Species, data = sf.alpha.div)##  significant p-value = 0.01
FSA::dunnTest(Observed~Species, data = sf.alpha.div) # signif between S.notata and S.porcus P value = 0.02

# Faith 
shapiro.test(sf.alpha.div$PD) # p-value = 1.379e-09
kruskal.test(PD~Species, data = sf.alpha.div)##  significant p-value = 0.03
FSA::dunnTest(PD~Species, data = sf.alpha.div) # Not signif 


# BOXPLOT

# highlight point from study Rodriguez- Ruiz in green

microfish<-sf.alpha.div %>% filter(FishID == "Microfish_3")
kormas<-sf.alpha.div%>% subset(FishID %in% c("SPO10.20" ,"SPO2.17", "SPO6.18","SPO9.19","SSC01.5","SSC02.6", "SSC04.7","SSC05.8"))

my_comparisons<-list(c("Snotata","Sporcus"),c("Snotata","Sscrofa"),c("Sporcus","Sscrofa"))
sf.alpha.div%>% group_by(Species)%>%summarise(freq=n())

ggp.shannon.Species<-
  sf.alpha.div %>% 
  ggplot(aes(Species,Shannon))+
  geom_boxplot(data=sf.alpha.div, fill= species_pal2, aes(alpha=0.5))+
  geom_point()+
  geom_point(data= microfish, aes(Species,Shannon), color ='green', alpha =0.3, size=4)+
  geom_point(data= kormas, aes(Species,Shannon), color ='red', alpha=0.3,size=4)+
  stat_compare_means(comparisons = my_comparisons,(aes(label = ..p.signif..)))+
  #scale_x_discrete(labels=c("1" = "Banyuls-sur-mer", "2" = "Carry-Le-Rouet", "3"="Bonifacio"))+
  ylab("Shannon's index")+ xlab("")+ylim(0,8)+
  # stat_summary(fun.y="mean", color ="red", shape= 20)+
  theme(legend.key.size = unit(4, "mm"), # resize legend so it fits in the figure
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45,hjust=0.9, size = 18, face = "italic"))

ggp.shannon.Species


ggarrange(ggp.shannon.SpeciesBO.glom,ggp.shannon.Species)



## 1.3. Beta diversity-----------------------------------------------------------

# First perform the analyses at the ASV level for the Species dataset 

# BO

# select the samples that are in the rarefied object
samples.BO<-species.meta.BO$FishID

#subset for the same individuals that are in the rarefied objects
BO.sub<-subset_samples(species.dataset.BO,FishID %in% samples.BO)

# run the function for the CLR transformation 
BO.clr.not.glom<-tran.clr.not.glom(BO.sub) 

# make the pca

pcaSpecies(phy.clr = BO.clr.not.glom, palette = species_pal2)

# Test variance 
permdisp(BO.clr.not.glom,"Species") # pvalue 0.122

#Test 

WDS.BO.not.glom<-Welch.Manova(BO.clr.not.glom,"Species") # p.value = 0.006 ; statistics = 1.29

res<-pairwise.adonis(BO.clr.not.glom,"Species")
# 
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 Snotata vs Sporcus  1  7204.014 1.284469 0.03640323   0.063      0.189    
# 2 Snotata vs Sscrofa  1 11270.131 1.449108 0.07086412   0.007      0.021   .
# 3 Sporcus vs Sscrofa  1 12236.848 1.849906 0.05308210   0.001      0.003   *


# Agglomerate the ASVs at the Genus level 

# Then I will agglomerate the phy obj from BO and from GY
# agglomerate because I want to compare the results with those I will obtain 
# comparing the three species across the Mediterranean; in this last case I will include 
# data from two other studies and so I will need to agglomerate to the genus level to be more conservative 

# BO
# select the samples that are in the rarefied object 
samples.BO<-species.meta.BO$FishID

#subset for the same individuals that are in the rarefied objects
BO.sub<-subset_samples(species.dataset.BO,FishID %in% samples.BO)

# run the function 
BO.clr<-tran.clr(BO.sub) 

# make the pca

pcaSpecies(phy.clr = BO.clr, palette = species_pal2)

# Test variance 
permdisp(BO.clr,"Species") # pvalue 0.104

#Test 

WDS.BO<-Welch.Manova(BO.clr,"Species") # p.value =  0.073 ; statistics =  1.25


# GY

# select the samples that are in the rarefied object 
samples.GY<-species.meta.GY$FishID

#subset for the same individuals that are in the rarefied objects
GY.sub<-subset_samples(species.dataset.GY,FishID %in% samples.GY)

# run the function 
GY.clr<-tran.clr(GY.sub) 

# make the pca

pcaSpecies(phy.clr = GY.clr, palette = species_pal2)

# Test variance 
permdisp(GY.clr,"Species") # pvalue 0.507

#Test 
permanova.GY<-permanova(GY.clr,"Species") # p.value=0.233

# 2 - Results chapter: "Phylosymbiosis In Mediterranean Scorpaena Fishes" ####
library(TreeDist)
library(ape)

#Load the matrix with genetic distances obtained from Turan et al., 2009 - 16S loci
gen.dist16<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_batch1_and_batch2_R1/mantel test/genetic_distances_Turan_applied_my_data.csv")
names<-gen.dist16$X
colnames(gen.dist16)<-c("",names)

rownames(gen.dist16)<-gen.dist16[,1]

gen.dist16<-gen.dist16[,-1]

gen.dist16.mat<-sapply(gen.dist16, as.numeric)

gen.dist16.mat<-as.matrix(gen.dist16) %>% as.dist()

# Generate the dissimilarity matrix - ASV level
clr.matrix.BO<-phyloseq::distance(BO.clr.not.glom,method="euclidean") #generate a distance matrix

#run mantel test
set.seed(1000)
mantel16<-vegan::mantel(xdis=clr.matrix.BO,ydis=gen.dist16.mat,method="pearson",permutations = 9999)

# Dissimilarity plot - pairwise scatter plot 
# from the tutorial https://jkzorz.github.io/2019/07/09/scatter-plots.html 

g16<-as.vector(gen.dist16.mat)
bac<-as.vector(clr.matrix.BO) 

new.data16<-data.frame(g16,bac)

#plot 
mm16 = ggplot(new.data16, aes(y = bac, x = g16,bac)) + 
  geom_point(size = 4, alpha = 0.1) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2)+  
  #scale_color_manual(values = color_by_ID)+
  labs(x = "Genetic distance - 16S", y = "Microbial dissimilarity (Aitchinson's distance)") + 
  scale_y_continuous(breaks = c(50,100,150))+
  scale_x_continuous(breaks=c(0.00,0.007,0.073,0.076))+
  ggtitle("r = 0.26 p.val = 0.0095", )+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10,angle=45, hjust=0.9), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"),
         plot.title = element_text(hjust = 0.5, size = 25),
         panel.background = element_blank())

mm16


#3 - Results chapter: "The effect of geographic location on the gut microbiota" ####

# Here data are agglomerated at the genus level, since they come from three different studies

# Subset the phyloseq object for the Species and Geo dataset
geo.dataset.sp<-subset_samples(sf, Species %in% "Sporcus")
geo.dataset.sn<-subset_samples(sf, Species %in% "Snotata")
geo.dataset.ss<-subset_samples(sf, Species %in% "Sscrofa")

# Rarefy 
geo.rare.sp<-tran.rare(geo.dataset.sp, sample.size=18025) # Keep 18025 because one of the external samples has a lower number of reads than 20000
geo.rare.sn<-tran.rare(geo.dataset.sn, sample.size=20000)
geo.rare.ss<-tran.rare(geo.dataset.ss, sample.size=18025) 

#Agglomerate to the genus level the three objects whit the external data 
geo.glom.sp<-tax_glom(geo.rare.sp, taxrank = "Genus")
geo.glom.sn<-tax_glom(geo.rare.sn, taxrank = "Genus")
geo.glom.ss<-tax_glom(geo.rare.ss, taxrank = "Genus")

#Obtain the metadata from each subset 
geo.sp.meta<-sample_data(geo.rare.sp) %>% data.frame()
geo.sn.meta<-sample_data(geo.rare.sn) %>% data.frame()
geo.ss.meta<-sample_data(geo.rare.ss) %>% data.frame()

## 3.1. Size differences among S.notata-----------------------------------------------

# Test the difference of sizes in the three locations
geo.sn.meta$MPA<-as.character(geo.sn.meta$MPA)
ggplot(geo.sn.meta,aes(MPA,Length, fill=MPA))+
  geom_boxplot()+
  geom_point()+
  stat_compare_means()+
  stat_compare_means()+
  scale_fill_manual(values = mpa.palette.sn)+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        legend.text = element_text(size=25))

# Test the difference of sizes in the three locations

shapiro.test(geo.sn.meta$Length) # not normally distributed
kruskal.test(geo.sn.meta$Length~geo.sn.meta$MPA) # p-value = 7.631e-06
FSA::dunnTest(Length~MPA, data=geo.sn.meta)

# Comparison          Z      P.unadj        P.adj
# 1      1 - 2  4.7237184 2.315711e-06 6.947132e-06
# 2      1 - 3  2.6498832 8.051960e-03 1.610392e-02
# 3      2 - 3 -0.6849235 4.933922e-01 4.933922e-01

## 3.2. Alpha diversity-----------------------------------------------------------------------

geo.sp.alpha.div<-alpha.div(geo.glom.sp,geo.sp.meta)
geo.sn.alpha.div<-alpha.div(geo.glom.sn,geo.sn.meta)
geo.ss.alpha.div<-alpha.div(geo.glom.ss,geo.ss.meta)


#S.notata
#Shannon
geo.sn.alpha.div$MPA<-as.factor(geo.sn.alpha.div$MPA)
shapiro.test(geo.sn.alpha.div$Shannon) # p-value = 0.002225
kruskal.test(Shannon~MPA, data = geo.sn.alpha.div) ## significant p-value = 0.001544

# Post Hoc -Tukey or Dunn
FSA::dunnTest(Shannon~MPA, data = geo.sn.alpha.div)

# 
# Comparison      Z     P.unadj      P.adj
#   1 - 2  -3.2700425 0.001075313 0.00322594
#   1 - 3  -2.5012251 0.012376447 0.02475289
#   2 - 3  -0.1923208 0.847490893 0.84749089


# Observed
shapiro.test(geo.sn.alpha.div$Observed) # p-value = 4.638e-08
kruskal.test(Observed~MPA, data = geo.sn.alpha.div) ## significant p-value = 4.419e-06

# Post Hoc -Tukey or Dunn
FSA::dunnTest(Observed~MPA, data = geo.sn.alpha.div)

# Comparison         Z      P.unadj        P.adj
# 1      1 - 2 -3.678284 2.348088e-04 4.696176e-04
# 2      1 - 3 -4.367720 1.255504e-05 3.766511e-05
# 3      2 - 3 -1.748318 8.040901e-02 8.040901e-02

# Faith 
shapiro.test(geo.sn.alpha.div$PD) # p-value = 5.405e-06
kruskal.test(PD~MPA, data = geo.sn.alpha.div) ## significant p-value = 4.633e-06

# Post Hoc -Tukey or Dunn
FSA::dunnTest(PD~MPA, data = geo.sn.alpha.div)

# Comparison         Z      P.unadj        P.adj
# 1      1 - 2 -3.723402 1.965561e-04 3.931123e-04 
# 2      1 - 3 -4.321488 1.549802e-05 4.649405e-05
# 3      2 - 3 -1.671385 9.464566e-02 9.464566e-02


#S.porcus 

#Shannon
geo.sp.alpha.div$MPA<-as.factor(geo.sp.alpha.div$MPA)
shapiro.test(geo.sp.alpha.div$Shannon) # p-value = 0.007756
wilcox.test(Shannon~MPA, data=geo.sp.alpha.div, paired = FALSE) # p-value = 0.341

#Observed 
shapiro.test(geo.sp.alpha.div$Observed) # p-value = 1.336e-05
wilcox.test(Observed~MPA, data=geo.sp.alpha.div, paired = FALSE) # p-value = 0.07651

#Faith 
shapiro.test(geo.sp.alpha.div$PD) # p-value = 0.0001661
wilcox.test(PD~MPA, data=geo.sp.alpha.div, paired = FALSE) # p-value = 0.05128


# S.scrofa

#Shannon
geo.ss.alpha.div$MPA<-as.factor(geo.ss.alpha.div$MPA)

wilcox.test(Shannon~MPA, data=geo.ss.alpha.div, paired = FALSE) # p-value = 0.3736

#Observed 
shapiro.test(geo.ss.alpha.div$Observed) # p-value = 0.006001
wilcox.test(Observed~MPA, data=geo.ss.alpha.div, paired = FALSE) # p-value = 0.1786

#Faith 
shapiro.test(geo.ss.alpha.div$PD) # p-value = 0.02387
wilcox.test(PD~MPA, data=geo.ss.alpha.div, paired = FALSE) # p-value = 0.1419



# PLOTS
# highlight point from study Rodriguez- Ruiz in green
# points

microfish<-geo.sn.alpha.div %>% filter(FishID == "Microfish_3")
kormas.sp<-geo.sp.alpha.div%>% subset(FishID %in% c("SPO10.20" ,"SPO2.17", "SPO6.18","SPO9.19"))
kormas.ss<-geo.ss.alpha.div%>% subset(FishID %in% c("SSC01.5","SSC02.6", "SSC04.7","SSC05.8"))



# S.notata 

geo.sn.plot<-alpha.plot(geo.sn.alpha.div,palette=mpa.palette.sn,points=microfish,color.points="green", method = "anova")
geo.sp.plot<-alpha.plot(geo.sp.alpha.div,palette=mpa.palette.sp.ss,points=kormas.sp,color.points="red",method = "wilcox.test")
geo.ss.plot<-alpha.plot(geo.ss.alpha.div,palette=mpa.palette.sp.ss,points=kormas.ss,color.points="red",method = "wilcox.test")


## 3.3 Beta diversity------------------------------------------------------------

# S.notata 

# select the samples that are in the rarefied object 
samples.sn<-geo.sn.meta$FishID

#subset for the same individuals that are in the rarefied objects
geo.sub.sn<-subset_samples(geo.dataset.sn,FishID %in% samples.sn)

# run the function - ASVs agglomerated at the Genus level 
geo.clr.sn<-tran.clr(geo.sub.sn)

# make the pca

pca.SN<-pcaMPA(phy.clr = geo.clr.sn, palette = mpa.palette.sn)

# Test variance 
set.seed(1000)
permdisp(geo.clr.sn,"MPA") # pvalue 0.002

#Test 

WDS.geo.sn<-Welch.Manova(geo.clr.sn, "MPA") # p.value = 0.001 ; statistics = 3.13

#Pairwise adonis 
pairwise.geo.sn<-pairwise.adonis(geo.clr.sn,"MPA")

# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#   1 3 vs 1 1  6022.102 7.522330 0.14322156   0.001      0.003   *
#   2 3 vs 2  1  4514.131 3.457642 0.07777385   0.001      0.003   *
#   3 1 vs 2  1  2352.616 2.751152 0.04001609   0.001      0.003   *
#   



# S.notata 

# S.notata not agglomerated at Genus level and not external samples 
# To confirm that the same results are also obtained at the ASV level 

#subset for the same individuals that are from externam studies
geo.sub.my.sn<-subset_samples(geo.sub.sn,!FishID %in% "Microfish_3")

# run the function 
geo.clr.sn.no.glom<-tran.clr.not.glom(geo.sub.my.sn) 

# make the pca

pcaMPA(phy.clr = geo.clr.sn.no.glom, palette = mpa.palette.sn)

# Test variance 
set.seed(1000)
permdisp(geo.clr.sn.no.glom, "MPA") # pvalue 0.003

#Test 

WDS.geo.sn.no.glom<-Welch.Manova(geo.clr.sn.no.glom,"MPA") # p.value = 0.001 ; statistics = 3.512187

#Pairwise adonis 
pairwise.geo.sn.no.glom<-pairwise.adonis(geo.clr.sn.no.glom, "MPA")
# 
#  pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#  3 vs 1  1 23650.813 7.438062 0.14460231   0.001      0.003   *
#  3 vs 2  1 23729.935 5.317128 0.11479831   0.001      0.003   *
#  1 vs 2  1  8165.362 2.615344 0.03867974   0.001      0.003   *
  
# S.porcus - ASVs agglomerated to Genus level 

# select the samples that are in the rarefied object 
samples.sp<-geo.sp.meta$FishID

#subset for the same individuals that are in the rarefied objects
geo.sub.sp<-subset_samples(geo.dataset.sp,FishID %in% samples.sp)

# run the function 
geo.clr.sp<-tran.clr(geo.sub.sp) 

# make the pca

pca.SP<-pcaMPA(phy.clr = geo.clr.sp, palette = mpa.palette.sp.ss)

# Test variance 
set.seed(1000)
permdisp(geo.clr.sp,"MPA") # pvalue 0.087 

#Test 

WDS.geo.sp<-Welch.Manova(geo.clr.sp,"MPA") # p.value = 0.001 ; statistics = 5.33

# S.scrofa - MPA agglomerted at Genus level

# select the samples that are in the rarefied object 
samples.ss<-geo.ss.meta$FishID

#subset for the same individuals that are in the rarefied objects
geo.sub.ss<-subset_samples(geo.dataset.ss,FishID %in% samples.ss)

# run the function 
geo.clr.ss<-tran.clr(geo.sub.ss) 

# make the pca

pca.SS<-pcaMPA(phy.clr = geo.clr.ss, palette = mpa.palette.sp.ss)

# Test variance 
permdisp(geo.clr.ss,"MPA") # pvalue 0.064

#Test 
WDS.geo.ss<-Welch.Manova(geo.clr.ss,"MPA") # p.value = 0.003 ; statistics = 3.21

# Arrange the plots 

ggarrange(pca.SN,pca.SP,pca.SS, nrow=1, common.legend = TRUE)



# 4 - Effect of fishing location--------------------------------------------------

## Agglometated to Genus level

# TEST the effect of Fishing Zones in each species in each location

# S.notata

#subset for the same individuals that are in the rarefied objects
geo.sub.sn.noMicrofish<-subset_samples(geo.sub.sn,!FishID %in% "Microfish_3")

# subset in three locations for SN 
snBA<-subset_samples(geo.sub.sn.noMicrofish, MPA %in% 1)
snCR<-subset_samples(geo.sub.sn.noMicrofish, MPA %in% 2)
snBO<-subset_samples(geo.sub.sn.noMicrofish, MPA %in% 3)

# run the function 
snBA.clr<-tran.clr(snBA) 
snCR.clr<-tran.clr(snCR) 
snBO.clr<-tran.clr(snBO) 

# The Wd Test will not work if there are groups with only one observation
df<-sample_data(geo.sub.sn.noMicrofish)%>% data.frame()
summary.zones<-df %>% group_by(Zone_by_location)%>% summarise(freq=n())

#Zone 7 in BA and Zone 5 in BO have only one observation, delete those from the analyses.
snBA.clr.noZ7<-subset_samples(snBA.clr, !Zone %in% "Z7")
snBO.clr.noZ5<-subset_samples(snBO.clr, !Zone %in% "Z5")

# BA

permdisp(snBA.clr.noZ7,"Zone") # p value = 0.055
Welch.snBA<-Welch.Manova(snBA.clr.noZ7,"Zone") # F value = 0.97 p value = 0.047

pairwise.snBA<-pairwise.adonis(snBA.clr.noZ7,"Zone") 

# pairs Df SumsOfSqs   F.Model         R2   p.value p.adjusted sig
# 1  Z1 vs Z11  1 1061.6084 1.4765860 0.19749469 0.1700000      1.000    
# 2   Z1 vs Z2  1 1021.2474 2.3137915 0.22433957 0.0580000      0.870    
# 3   Z1 vs Z4  1 1057.4249 3.8705917 0.32606561 0.0230000      0.345    
# 4   Z1 vs Z6  1 1113.7867 3.3281154 0.29379251 0.0430000      0.645    
# 5   Z1 vs Z8  1 1233.5581 1.3315392 0.39967688 0.3333333      1.000    
# 6  Z11 vs Z2  1  765.9716 1.5242116 0.11270244 0.0380000      0.570    
# 7  Z11 vs Z4  1 1153.8527 2.9554389 0.19761633 0.0020000      0.030 *****
# 8  Z11 vs Z6  1  931.0896 2.1583299 0.15244241 0.0050000      0.075    
# 9  Z11 vs Z8  1 1134.3782 1.5638422 0.20675235 0.1090000      1.000    
# 10  Z2 vs Z4  1  276.2159 0.9909793 0.06610504 0.4340000      1.000    
# 11  Z2 vs Z6  1  224.4776 0.7152311 0.04860482 0.8470000      1.000    
# 12  Z2 vs Z8  1 1405.1217 3.1491759 0.28245818 0.0190000      0.285    
# 13  Z4 vs Z6  1  299.2579 1.3743137 0.08939025 0.1140000      1.000    
# 14  Z4 vs Z8  1 1610.4885 5.7929620 0.41999405 0.0200000      0.300    
# 15  Z6 vs Z8  1 1508.7620 4.4444241 0.35714180 0.0280000      0.420  
# 

# CR
permdisp(snCR.clr,"Zone") # p value = 0.17
permanova.snCR<-permanova(snCR.clr,"Zone") # F value = 1.11 p value =0.19

# BO
permdisp(snBO.clr.noZ5,"Zone") # p value = 0.123

permanova.snBO<-permanova(snBO.clr.noZ5,"Zone") # F value = 1.11 p value =0.19

pairwise.snBO<-pairwise.adonis(snBO.clr.noZ5,"Zone") # F value = 0.98 p value =0.049 

#  pairs Df SumsOfSqs  F.Model        R2   p.value p.adjusted sig
# 1 Z1 vs Z2_3  1  1984.596 1.184444 0.1648624 0.2080000      0.624    
# 2   Z1 vs Z6  1  1625.743 1.033176 0.1469004 0.2840000      0.852    
# 3 Z2_3 vs Z6  1  1393.260 1.270917 0.3885508 0.3333333      1.000

pca.snBO<-pcaZone(phy.clr = snBO.clr.noZ5, palette = zone.palette)

# make the pca
set.seed(1000)
pca.snBA<-pcaZone(phy.clr = snBA.clr.noZ7, palette = zone.palette)
pca.snCR<-pcaZone(phy.clr = snCR.clr, palette = zone.palette)
pca.snBO<-pcaZone(phy.clr = snBO.clr.noZ5, palette = zone.palette)


# S. porcus

spBO.sub<-subset_samples(BO.sub,Species %in% "Sporcus")
spBO.clr<-tran.clr(spBO.sub)

df.spBO<-sample_data(spBO.clr)%>% data.frame()
summary.zones.sp<-df.spBO %>% group_by(Zone_by_location)%>% summarise(freq=n())
spBO.clr.noZ5Z6<-subset_samples(spBO.clr, !Zone %in% c("Z5","Z6"))

permdisp(spBO.clr.noZ5Z6,"Zone") # p value = 0.032
WDS.location.spBO<-Welch.Manova(spBO.clr.noZ5Z6,"Zone") # F value = 1.44  p value = 0.001
pairwise.spBO<-pairwise.adonis(spBO.clr.noZ5Z6,"Zone") # F value = 0.98 p value =0.049 

pca.spBO<-pcaZone(phy.clr = spBO.clr.noZ5Z6, palette = zone.palette)

# S. scrofa

ssBO.sub<-subset_samples(BO.sub,Species %in% "Sscrofa")
ssBO.clr<-tran.clr(ssBO.sub)
WDS.location.ssBO<-Welch.Manova(ssBO.clr,"Zone") # F value = 0.798  p value = 0.515


# Not agglomerated 
# run the function 
# S.notata
snBA.clr.no.glom<-tran.clr.not.glom(snBA) 
snCR.clr.no.glom<-tran.clr.not.glom(snCR) 
snBO.clr.no.glom<-tran.clr.not.glom(snBO) 

sample_data(snBO.clr.no.glom)%>% data.frame()%>%group_by(Zone)%>%summarise(freq=n())

#Zone 7 in BA and Zone 5 in BO have only one observation, delete those from the analyses.
snBA.clr.no.glom.noZ7<-subset_samples(snBA.clr.no.glom, !Zone %in% "Z7")
snBO.clr.no.glom.noZ5<-subset_samples(snBO.clr.no.glom, !Zone %in% "Z5")

# BA
# Test
permdisp(snBA.clr.no.glom.noZ7,"Zone") # p value = 0.052
Welch.snBA.no.glom<-Welch.Manova(snBA.clr.no.glom.noZ7,"Zone") # F value = 1.01 p value = 0.008
pairwise.snBA.no.glom<-pairwise.adonis(snBA.clr.no.glom.noZ7,"Zone") 

# pairs Df SumsOfSqs   F.Model         R2   p.value p.adjusted sig
# 1  Z1 vs Z11  1  3645.373 1.2727975 0.17500796 0.2110000      1.000    
# 2   Z1 vs Z2  1  4142.580 2.1769564 0.21391036 0.0140000      0.210    
# 3   Z1 vs Z4  1  4084.719 2.7234063 0.25396840 0.0230000      0.345    
# 4   Z1 vs Z6  1  3947.655 2.2651694 0.22066556 0.0200000      0.300    
# 5   Z1 vs Z8  1  3858.178 0.9560997 0.32343283 0.6666667      1.000    
# 6  Z11 vs Z2  1  5619.420 2.6634531 0.18163887 0.0020000      0.030 ***
# 7  Z11 vs Z4  1  5988.647 3.2527298 0.21325558 0.0010000      0.015 ***
# 8  Z11 vs Z6  1  5161.991 2.5770626 0.17678888 0.0010000      0.015 ***
# 9  Z11 vs Z8  1  3907.965 1.2908098 0.17704615 0.1800000      1.000    
# 10  Z2 vs Z4  1  1735.959 1.2071814 0.07938232 0.1350000      1.000    
# 11  Z2 vs Z6  1  1619.720 1.0271995 0.06835602 0.3460000      1.000    
# 12  Z2 vs Z8  1  5987.046 2.9558036 0.26979341 0.0190000      0.285    
# 13  Z4 vs Z6  1  1396.418 1.0370674 0.06896739 0.3570000      1.000    
# 14  Z4 vs Z8  1  6077.578 3.7459119 0.31891197 0.0200000      0.300    
# 15  Z6 vs Z8  1  5512.394 2.9551319 0.26974864 0.0280000      0.420    


#CR
permdisp(snCR.clr.no.glom,"Zone") # p value = 0.232
permanova.snCR.no.glom<-permanova(snCR.clr.no.glom,"Zone") # F value = 1.13 p value =0.352

# BO
permdisp(snBO.clr.no.glom.noZ5,"Zone") # p value =0.146
permanova.snBO.no.glom<-permanova(snBO.clr.no.glom.noZ5,"Zone") #  p value = 0.332

pairwise.snBO<-pairwise.adonis(snBO.clr.no.glom.noZ5,"Zone") # F value = 0.98 p value =0.049 

# pairs Df SumsOfSqs  F.Model        R2   p.value p.adjusted sig
# 1 Z1 vs Z2_3  1  6948.202 1.110607 0.1561902 0.3480000          1    
# 2   Z1 vs Z6  1  6637.301 1.120830 0.1574016 0.4650000          1    
# 3 Z2_3 vs Z6  1  5139.714 1.175257 0.3701297 0.3333333          1    



# make the pca

pca.snBA<-pcaZone(phy.clr = snBA.clr.noZ7, palette = zone.palette)
pca.snCR<-pcaZone(phy.clr = snCR.clr, palette = zone.palette)
pca.snBO<-pcaZone(phy.clr = snBO.clr.noZ5, palette = zone.palette)


# S. porcus

spBO.clr.no.glom<-tran.clr.not.glom(spBO.sub)

spBO.clr.not.glom.noZ5Z6<-subset_samples(spBO.clr.no.glom, !Zone %in% c("Z5","Z6"))

permdisp(spBO.clr.not.glom.noZ5Z6,"Zone") # p value = 0.029 
WDS.location.spBO.not.glom<-Welch.Manova(spBO.clr.not.glom.noZ5Z6,"Zone") # F value = 1.19 p value = 0.003

pca.spBO.not.glom<-pcaZone(phy.clr = spBO.clr.not.glom.noZ5Z6, palette = zone.palette)
pairwise.spBO.not.glom<-pairwise.adonis(spBO.clr.not.glom.noZ5Z6,"Zone") # F value = 0.98 p value =0.049 

# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1   Z1 vs Z11  1  9343.308 1.3757191 0.12093470   0.017      0.102    
# 2   Z1 vs Z13  1  4108.571 0.9601378 0.12061824   0.280      1.000    
# 3  Z1 vs Z2_3  1  5750.165 1.3390308 0.08729566   0.045      0.270    
# 4  Z11 vs Z13  1  6446.524 0.7724072 0.13381025   0.768      1.000    
# 5 Z11 vs Z2_3  1 11292.542 1.8848158 0.13574655   0.008      0.048   .
# 6 Z13 vs Z2_3  1  4210.905 1.1168160 0.11039204   0.331      1.000 

# S. scrofa
ssBO.clr.not.glom<-tran.clr.not.glom(ssBO.sub)
WDS.location.ssBO.not.glom<-Welch.Manova(ssBO.clr.not.glom,"Zone") # F value = 0.862  p value = 0.165


# 5. Effect of the size of the fishes------------------------------------------------------------------------------

## PERMANOVA for the size - NOT agglomerated to genus level
# S.notata 
set.seed(1000)
permanova.size.SN.BO<-permanova(snBO.clr.no.glom,"Length") # stat = 0.084, P value = 0.687

# S.porucs 
set.seed(1000)
permanova.size.SP.not.glom<-permanova(spBO.clr.no.glom,"Length") # stat = 0.732, P value = 0.905

# S.scrofa 
set.seed(1000)
permanova.size.SS.not.glom<-permanova(ssBO.clr.not.glom,"Length") # stat = 0.07259, P value = 0.893



## PERMANOVA for the size - agglomerated to genus level

### S.notata 
set.seed(1000)
geo.clr.sn.noMicrofish<-subset_samples(geo.clr.sn,!FishID %in% "Microfish_3")
permanova.size.SN<-permanova(geo.clr.sn.noMicrofish,"Length") # stat = 0.95427, P value = 0.585

### S.porucs 
set.seed(1000)
permanova.size.SP<-permanova(spBO.clr,"Length") #  P value = 0.88

### S.scrofa 
set.seed(1000)
permanova.size.SS<-permanova(ssBO.clr,"Length") #  P value = 0.874


# 6. Alpha diversity using the Extended dataset

sf

# Rarefy 
sf.rare<-tran.rare(sf, sample.size=18025) # Keep 18025 because one of the external samples has a lower number of reads than 20000

#Agglomerate to the genus level the three objects whit the external data 
sf.rare.glom<-tax_glom(sf.rare, taxrank = "Genus", NArm = FALSE)

#Obtain the metadata from each subset 
sf.meta<-sample_data(sf.rare) %>% data.frame()


# Genrate the alpha divertiy table 
sf.alpha.div<-alpha.div(sf.rare.glom,sf.meta)

# Shannon 
shapiro.test(sf.alpha.div$Shannon) # p-value = 5.624e-08
kruskal.test(Shannon~Species, data = sf.alpha.div)##  significant p-value = 0.007
FSA::dunnTest(Shannon~Species, data = sf.alpha.div) ## siginf between S.scrofa and S.notata p value = 0.013

# Obsrved 
shapiro.test(sf.alpha.div$Observed) # p-value =9.044e-12
kruskal.test(Observed~Species, data = sf.alpha.div)##  significant p-value = 0.01
FSA::dunnTest(Observed~Species, data = sf.alpha.div) # signif between S.notata and S.porcus P value = 0.02

# Faith 
shapiro.test(sf.alpha.div$PD) # p-value = 1.379e-09
kruskal.test(PD~Species, data = sf.alpha.div)##  significant p-value = 0.03
FSA::dunnTest(PD~Species, data = sf.alpha.div) # Not signif

# Boxplot

# highlight point from study Rodriguez- Ruiz in green

microfish<-sf.alpha.div %>% filter(FishID == "Microfish_3")
kormas<-sf.alpha.div%>% subset(FishID %in% c("SPO10.20" ,"SPO2.17", "SPO6.18","SPO9.19","SSC01.5","SSC02.6", "SSC04.7","SSC05.8"))

my_comparisons<-list(c("Snotata","Sporcus"),c("Snotata","Sscrofa"),c("Sporcus","Sscrofa"))

ggp.shannon.Species<-
  sf.alpha.div %>% 
  ggplot(aes(Species,Shannon))+
  geom_boxplot(data=sf.alpha.div, fill= species_pal2, aes(alpha=0.5))+
  geom_point()+
  geom_point(data= microfish, aes(Species,Shannon), color ='green', alpha =0.3, size=4)+
  geom_point(data= kormas, aes(Species,Shannon), color ='red', alpha=0.3,size=4)+
  stat_compare_means(comparisons = my_comparisons,(aes(label = ..p.signif..)))+
  #scale_x_discrete(labels=c("1" = "Banyuls-sur-mer", "2" = "Carry-Le-Rouet", "3"="Bonifacio"))+
  ylab("Shannon's index")+ xlab("")+ylim(0,8)+
  # stat_summary(fun.y="mean", color ="red", shape= 20)+
  theme(legend.key.size = unit(4, "mm"), # resize legend so it fits in the figure
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45,hjust=0.9, size = 18, face = "italic"))

ggp.shannon.Species


## NOTE 
# Eventually we cannot trust these results because S.notata displays significant different values of alpha diversity indexes among the three locations in which it was collected (see point 3.2)



# 7.- Barplot Top genera - Fig 1C----------------------------------------------------

#Create functions:

#Function to agglomerate the ASVs at a taxonomic level
taxglom_topN<-function(ps.object,taxrank,n=topX){
  top<-tax_glom(ps.object,taxrank=taxrank, NArm = FALSE)
  top_names<-names(sort(taxa_sums(top),decreasing = TRUE)[1:n])
  top_glom<-prune_taxa(top_names,top)
}


#Create a function to build stacked barplots at different taxonomical levels with the top nASVs
stacked.barplot.topASVs<-function(ps.object,taxrank,x.axis.value, grid,n=topX, palette){
  ps.object.glom<-tax_glom(ps.object,taxrank = taxrank)
  topX<-names(sort(taxa_sums(ps.object.glom), decreasing=TRUE))[1:n]
  topX_otus <- prune_taxa(topX, ps.object.glom)
  Ntop<-transform_sample_counts(topX_otus,function(x) {x/sum(x)})
  p=plot_bar(Ntop,x=x.axis.value,fill=taxrank)
  p=p+geom_bar(position="fill",stat = "identity")
  p=p+scale_fill_manual(values = palette)
  p=p+facet_wrap(grid,scale="free_x")
  p=p+theme(legend.key.size = unit(3, "mm"), 
            legend.position = c("bottom"), 
            legend.box = c("horizontal"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(p)
}

# Dataset used: Species dataset rarefied
species.rare.BO

# Calculate the relative abundances of each taxa in the otu table 

species.rare.BO.ra<-transform_sample_counts(species.rare.BO,function(x) {x/sum(x)})
speciesBO.glom<-taxglom_topN(species.rare.BO.ra, taxrank = "Genus", n=12)
speciesBO.otu<-as.data.frame(t(otu_table(speciesBO.glom)))
speciesBO.tax<-as.data.frame(tax_table(speciesBO.glom))
top10.sum<-colSums(speciesBO.otu[,c(1:46)]) #relative abundance of the Gamma families top 10
other<-1-top10.sum[1:46] #obtain the relative abundance of all the other families

#import the other row in the otu table
speciesBO.otu.otu.ed<-rbind(speciesBO.otu,other)
rownames(speciesBO.otu.otu.ed)[13]<-"Zother genera"
speciesBO.otu.otu.ed<-data.matrix(t(speciesBO.otu.otu.ed))

#import the other row in the taxa table
speciesBO.tax.ed<-speciesBO.tax
speciesBO.tax.ed[13,]<-c("Zother genera","Zother genera","Zother genera","Zother genera","Zother genera","Zother genera") #"Zother genera", remove the Z in post processing of the image. I did this just to positionate the other genera at the bottom of the bars
rownames(speciesBO.tax.ed)[13]<-"Zother genera"
taxa_names<-rownames(speciesBO.tax.ed)
speciesBO.tax.ed<-as.matrix(speciesBO.tax.ed)

#recombine the phylo obj

speciesBO.glom.ed<- phyloseq(otu_table(speciesBO.otu.otu.ed, taxa_are_rows = FALSE),
                             tax_table(speciesBO.tax.ed), taxa_names(taxa_names),
                             sample_data(sample_data(speciesBO.glom)))

install.packages("ggsci")  
library(RColorBrewer)
library(ggsci)


bar12Genera<-plot_bar(speciesBO.glom.ed,fill="Genus",x="WP3_alldata_FishID_Nanodrop")+
  geom_bar(position="fill",stat = "identity")+
  facet_wrap(~Species,scale="free_x")+
  scale_fill_igv()+
  xlab("Sample ID")+ ylab("Relative abundance")+
  theme(legend.key.size = unit(3, "mm"), 
        legend.position = c("bottom"), 
        legend.box = c("horizontal"),
        legend.text = element_text(face="italic", size= 20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# 8 - Differential abundance ----------------------------------------------------
## 8.1 Geographic location-----------------------------------------------------------

# S. notata 

library(DESeq2)

### Set up the DESeq function------------------------------------------------------------------------

# Normally DESeq normalize raw counts by log-based scaling automatically. 
# but
# DESeq cannot compute log geometric means when every gene contains zeros
# Therefore
# work around in https://github.com/joey711/phyloseq/issues/387
# calculate size factors using a zero-tolerant variant of geometric mean
# included in package vignette

# Then

# Taxa-wise variance is estimated. These values tell how much each taxa varies between samples.
# A curve is fitted over all those taxa-wise variance estimates that we got in the last step.
# This model tells how big the variance is in a specific abundance level.
# The model is used to shrink those individual variance estimates to avoid the effect of, e.g., small sample size and higher variance. 
# This reduces the likelihood to get false positives.
# Variance estimates are used to compare different groups. We receive a result that shows whether the variance is explained by groups



gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


deseq.fun.2<-function(pair, ref= c()){
  ds2<-phyloseq_to_deseq2(pair,~MPA)
  
  geoMeans = apply(counts(ds2), 1, gm_mean)
  ds2 = estimateSizeFactors(ds2, geoMeans = geoMeans)
  
  ds2$Species<-relevel(ds2$MPA, ref = ref) # in this way you reset your reference variable, so the abundances will be compared to the furthest zone
  
  # Run the DeSeq fun
  dds <- DESeq(ds2, fitType="local")
  
  # Gets the results from the object
  res <- results(dds)
  
  # Creates a data frame from results
  df <- as.data.frame(res)
  
  # Adds taxon column that includes names of taxa
  df <- rownames_to_column(df,var="ASV")
  t<-data.frame(tax_table(pair))%>% rownames_to_column(var="ASV")
  df<-left_join(df,t,by="ASV")
  
  # Orders the rows of data frame in increasing order firstly based on column
  # "log2FoldChange" and secondly based on "padj" column
  df <- df %>% arrange(log2FoldChange, padj)
  
  #kable(head(df)) %>% 
  # kable_styling("striped") %>% 
  #scroll_box(width = "100%")
  
}

# select the samples that are in the rarefied object 
samples.sn<-geo.sn.meta$FishID

#subset for the same individuals that are in the rarefied objects
geo.sub.sn<-subset_samples(geo.dataset.sn,FishID %in% samples.sn)

# Keep only the S. notata from this study
geo.dataset.sn.Lilli<-subset_samples(geo.sub.sn, !FishID %in% "Microfish_3")

# remove the ASVs that have less than 20 counts 
geo.dataset.sn.Lilli<-prune_taxa(taxa_sums(geo.dataset.sn.Lilli) > 20, geo.dataset.sn.Lilli)

#check if the library size of the three groups is very unbalanced
counts1<-geo.dataset.sn.Lilli %>% subset_samples(MPA == 1) %>% sample_sums()
counts1<-as.data.frame(counts1)
counts1$MPA<-"MPA1"
colnames(counts1)<-c("counts","MPA")

counts2<-geo.dataset.sn.Lilli %>% subset_samples(MPA == 2) %>% sample_sums()
counts2<-as.data.frame(counts2)
counts2$MPA<-"MPA2"
colnames(counts2)<-c("counts","MPA")

counts3<-geo.dataset.sn.Lilli %>% subset_samples(MPA == 3) %>% sample_sums()
counts3<-as.data.frame(counts3)
counts3$MPA<-"MPA3"
colnames(counts3)<-c("counts","MPA")

counts.df<-rbind(counts1,counts2,counts3)

kruskal.test(counts~MPA, data=counts.df) #signif
boxplot(counts~MPA, data=counts.df)

# there is unequal library size and as suggested in Weiss et al 2017 I will rarefy the data before runnin Deseq2

sn.Lilli.rare<-tran.rare(geo.dataset.sn.Lilli, sample.size = 20000)

# run the function on the Genus level
sn.Lilli.rare.gen<-tax_glom(sn.Lilli.rare,taxrank = "Genus", NArm = FALSE) #661 taxa and 77 samples


# subset by couples of Species to run the fucntion on pairs
BA.CR.gen<-sn.Lilli.rare.gen%>% subset_samples(MPA %in% c(1,2)) #pair 1
CR.BO.gen<-sn.Lilli.rare.gen%>% subset_samples(MPA %in% c(2,3))  #pair 2                        
BA.BO.gen<-sn.Lilli.rare.gen%>% subset_samples(MPA %in% c(1,3))  #pair 3                        
sample_data(BA.CR.gen)$MPA<-as.factor(sample_data(BA.CR.gen)$MPA)
sample_data(CR.BO.gen)$MPA<-as.factor(sample_data(CR.BO.gen)$MPA)
sample_data(BA.BO.gen)$MPA<-as.factor(sample_data(BA.BO.gen)$MPA)

# Run DESeq2 on every pair
# specify the reference fish Species for each pair. This means that the logFoldChange2 > 0 will
# entail that the reference Species has a lower abundance of that taxa. LogFoldChange < 0 means that the reference species has higher abundance of the taxa.

p1.gen<-deseq.fun.2(BA.CR.gen,ref = c("1"))
p2.gen<-deseq.fun.2(CR.BO.gen,ref = c("2"))
p3.gen<-deseq.fun.2(BA.BO.gen,ref = c("1")) 

# Keep only the taxa that have padj < 0.05

p1.gen.flt<-subset(p1.gen, padj < 0.05)
p2.gen.flt<-subset(p2.gen, padj < 0.05)
p3.gen.flt<-subset(p3.gen, padj < 0.05)

### Set up the functions for ANCOM II----------------------------------------------------------------------------
library(exactRankTests)
library(nlme)
library(dplyr)

# OTU table should be a matrix/data.frame with each feature in rows and sample in columns. 
# Metadata should be a matrix/data.frame containing the sample identifier. 

# Data Pre-Processing 
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL, 
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1) 
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = residuals(f_fit)
    y = t(t(z) - e)
    
    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n
        
        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break
        
        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))
        
        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 + 
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }
      
      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = t(apply(y, 1, function(i) unlist(tapply(i, group, function(j) outlier_check(j)))))
    feature_table[out_ind] = NA
  }
  
  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }
  
  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }
  
  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    
    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
    
    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1
    
    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)
    
    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }
  
  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}

# ANCOM main function
ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH", 
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL){
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]
  
  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = exactRankTests::wilcox.exact
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }
  
  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted 
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ]) 
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata
    
    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        tfun(tformula, data = data.frame(x, alr_data, check.names = FALSE))$p.value
      }
      ) 
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE), 
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit)
        anova(fit)[main_var, "p-value"]
      }
      ) 
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  out = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out = out%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                     detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                     detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                     detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))
  
  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    res = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE, 
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE, 
                     row.names = NULL, check.names = FALSE)
    res[match(taxa_id, res$taxa_id), ] = out
  }else{
    res = out
  }
  
  return(res)
}


### RUN ANCOM II - for Geo location------------------------------------------------ 

#for this method I need to have the otu table with absolute values and not relative, not rarefied
#the table need to have the Samples as row and the taxa as col. The two tables need to have "Sample.ID" as the identifier of the samples.

ancom.fun<-function(pair){
  asv.tbl.ancom<-as.data.frame(t(otu_table(pair)))
  asv.tbl.ancom$otu_id<-rownames(asv.tbl.ancom)
  asv.tbl.ancom<-data.matrix(asv.tbl.ancom)
  metadata.ancom<-sample_data(pair) %>% data.frame()
  metadata.ancom$Sample.ID <-rownames(metadata.ancom)
  
  # Step 1: Data preprocessing
  
  feature_table <- asv.tbl.ancom
  meta_data<-metadata.ancom
  sample_var <- "Sample.ID"
  group_var <- NULL
  out_cut <- 0.05
  zero_cut <- 0.90
  lib_cut <- 1000
  neg_lb <- FALSE
  prepro <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                      out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # Step 2: ANCOM
  
  main_var<-"MPA"
  p_adj_method<-"BH"
  alpha<-0.05
  adj_formula<-NULL
  rand_formula<-NULL
  t_start<-Sys.time()
  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
              alpha, adj_formula, rand_formula)
  colnames(res)[1] <- "ASV"
  tax<-as.data.frame(tax_table(pair))    
  tax<- rownames_to_column(tax, var="ASV")
  res<-left_join(res,tax, by= "ASV")
  
}

##Apply to the different pair sets-GENERA

ban.carry.ancom<-ancom.fun(pair=BA.CR.gen)
ban.bon.ancom<-ancom.fun(pair=BA.BO.gen )
carry.bon.ancom<-ancom.fun(pair=CR.BO.gen)

ban.carry.ancom.sub<-subset(ban.carry.ancom,detected_0.7 == TRUE)
ban.bon.ancom.sub<-subset(ban.bon.ancom,detected_0.7 == TRUE)
carry.bon.ancom.sub<-subset(carry.bon.ancom,detected_0.7 == TRUE)


### Combine DESeq and ANCOM II (0.7) tables-------------------------------------
library(grid)
library(forcats)
#### BA-CR------------------------------------------------------------------------
dim(p1.gen.flt) # 6
dim(ban.carry.ancom.sub) # 5

BA.CR.combo<-inner_join(ban.carry.ancom.sub,p1.gen.flt,by="ASV")

# Create a plot to show the difference between these three Species at the Genus level

theme_set(theme_minimal())
# MPA 1 vs MPA 2 - reference : MPA 1
BA.CR.combo<-arrange(BA.CR.combo,log2FoldChange)
ba.cr.sn<-ggplot(BA.CR.combo, aes (x=fct_inorder(Genus.x), fill=log2FoldChange))+
  geom_bar()+
  geom_vline(xintercept=0, linetype='dotted', col = 'black', size =1)+
  scale_fill_gradient2(low="black", mid = "white",high="red" , midpoint = 0)+
  xlab("")+ ylab("")+
  theme(axis.text.x=element_text(angle = 45, hjust=1, size = 20, face= "italic"))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())


# annotate("segment", x=7, xend=6, y=1.2,
#### CR-BO-----------------------------------------------------------------------
dim(p2.gen.flt) # 13
dim(carry.bon.ancom.sub) # 5

CR.BO.combo<-inner_join(carry.bon.ancom.sub,p2.gen.flt,by="ASV")

cr.bo.sn<-ggplot(CR.BO.combo, aes (x=fct_inorder(Genus.x), fill=log2FoldChange))+
  geom_bar()+
  geom_vline(xintercept=0, linetype='dotted', col = 'black', size =1)+
  scale_fill_gradient2(low="black", mid = "white",high="red" , midpoint = 0)+
  xlab("")+ ylab("")+
  theme(axis.text.x=element_text(angle = 45, hjust=1, size = 20, face= "italic"))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())


#### BA-BO-----------------------------------------------------------------------
dim(p3.gen.flt) # 20
dim(ban.bon.ancom.sub) # 10

BA.BO.combo<-inner_join(ban.bon.ancom.sub,p3.gen.flt,by="ASV")

BA.BO.combo<-arrange(BA.BO.combo,log2FoldChange)
ba.bo.sn<-ggplot(BA.BO.combo, aes (x=fct_inorder(Genus.x), fill=log2FoldChange))+
  geom_bar()+
  geom_vline(xintercept=0, linetype='dotted', col = 'black', size =1)+
  scale_fill_gradient2(low="black", mid = "white",high="red" , midpoint = 0)+
  xlab("")+ ylab("")+
  theme(axis.text.x=element_text(angle = 45, hjust=1, size = 20, face= "italic"))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())

## 8.2. Species-----------------------------------------------------------------
deseq.fun.2.Species<-function(pair, ref= c()){
  ds2<-phyloseq_to_deseq2(pair,~Species)
  
  geoMeans = apply(counts(ds2), 1, gm_mean)
  ds2 = estimateSizeFactors(ds2, geoMeans = geoMeans)
  
  ds2$Species<-relevel(ds2$Species, ref = ref) # in this way you reset your reference variable, so the abundances will be compared to the furthest zone
  
  # Run the DeSeq fun
  dds <- DESeq(ds2, fitType="local")
  
  # Gets the results from the object
  res <- results(dds)
  
  # Creates a data frame from results
  df <- as.data.frame(res)
  
  # Adds taxon column that includes names of taxa
  df <- rownames_to_column(df,var="ASV")
  t<-data.frame(tax_table(pair))%>% rownames_to_column(var="ASV")
  df<-left_join(df,t,by="ASV")
  
  # Orders the rows of data frame in increasing order firstly based on column
  # "log2FoldChange" and secondly based on "padj" column
  df <- df %>% arrange(log2FoldChange, padj)
  
  #kable(head(df)) %>% 
  # kable_styling("striped") %>% 
  #scroll_box(width = "100%")
  
}

# Use the BO.sub object (only samples from BO from this study and already filtered for those with less than 20.000 reads)

BO.sub # unrarefied object 

#check if the library size of the three groups is very unbalanced - if it is as said in Weiss et al I should rarefy the data before running DESeq2
counts1.BO<-BO.sub %>% subset_samples(Species == "Snotata") %>% sample_sums()
counts1.BO<-as.data.frame(counts1.BO)
counts1.BO$Species<-"Snotata"
colnames(counts1.BO)<-c("counts","Species")

counts2.BO<-BO.sub %>% subset_samples(Species == "Sporcus") %>% sample_sums()
counts2.BO<-as.data.frame(counts2.BO)
counts2.BO$Species<-"Sporcus"
colnames(counts2.BO)<-c("counts","Species")

counts3.BO<-BO.sub %>% subset_samples(Species == "Sscrofa") %>% sample_sums()
counts3.BO<-as.data.frame(counts3.BO)
counts3.BO$Species<-"Sscrofa"
colnames(counts3.BO)<-c("counts","Species")

counts.df.BO<-rbind(counts1.BO,counts2.BO,counts3.BO)

kruskal.test(counts~Species, data=counts.df.BO) # not signif
# No need to rarefy before running DesEQ2

# run the function on the Genus level 
BO.sub.noZero<-nozero_taxa(BO.sub)
BO.sub.unrare.gen<-tax_glom(BO.sub.noZero,taxrank = "Genus", NArm = FALSE) 

# subset by couples of Species to run the fucntion on pairs
sn.ss.gen<-BO.sub.unrare.gen%>% subset_samples(Species %in% c("Snotata","Sscrofa")) #pair 1
sp.ss.gen<-BO.sub.unrare.gen%>% subset_samples(Species %in% c("Sporcus","Sscrofa"))  #pair 2                        
sp.sn.gen<-BO.sub.unrare.gen%>% subset_samples(Species %in% c("Sporcus","Snotata"))  #pair 3                        

# Run DESeq2 on every pair
# specify the reference fish Species for each pair. This means that the logFoldChange2 > 0 will
# entail that the reference Species has a lower abundance of that taxa. LogFoldChange < 0 means that the reference species has higher abundance of the taxa.

sn.ss<-deseq.fun.2.Species(sn.ss.gen,ref = c("Snotata"))
sp.ss<-deseq.fun.2.Species(sp.ss.gen,ref = c("Sscrofa"))
sp.sn<-deseq.fun.2.Species(sp.sn.gen,ref = c("Snotata")) 

# Keep only the taxa that have padj < 0.05

sn.ss.flt<-subset(sn.ss, padj < 0.05)
sp.ss.flt<-subset(sp.ss, padj < 0.05)
sp.sn.flt<-subset(sp.sn, padj < 0.05)

### RUN ANCOM II - for Species location------------------------------------------------ 

#for this method I need to have the otu table with absolute values and not relative, not rarefied
#the table need to have the Samples as row and the taxa as col. The two tables need to have "Sample.ID" as the identifier of the samples.

ancom.fun<-function(pair){
  asv.tbl.ancom<-as.data.frame(t(otu_table(pair)))
  asv.tbl.ancom$otu_id<-rownames(asv.tbl.ancom)
  asv.tbl.ancom<-data.matrix(asv.tbl.ancom)
  metadata.ancom<-sample_data(pair) %>% data.frame()
  metadata.ancom$Sample.ID <-rownames(metadata.ancom)
  
  # Step 1: Data preprocessing
  
  feature_table <- asv.tbl.ancom
  meta_data<-metadata.ancom
  sample_var <- "Sample.ID"
  group_var <- NULL
  out_cut <- 0.05
  zero_cut <- 0.90
  lib_cut <- 1000
  neg_lb <- FALSE
  prepro <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                      out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # Step 2: ANCOM
  
  main_var<-"Species"
  p_adj_method<-"BH"
  alpha<-0.05
  adj_formula<-NULL
  rand_formula<-NULL
  t_start<-Sys.time()
  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
              alpha, adj_formula, rand_formula)
  colnames(res)[1] <- "ASV"
  tax<-as.data.frame(tax_table(pair))    
  tax<- rownames_to_column(tax, var="ASV")
  res<-left_join(res,tax, by= "ASV")
  
}


##Apply to the different pair sets-GENERA

sn.ss.ancom<-ancom.fun(pair=sn.ss.gen)
sp.ss.ancom<-ancom.fun(pair=sp.ss.gen)
sp.sn.ancom<-ancom.fun(pair=sp.sn.gen)

sn.ss.ancom.sub<-subset(sn.ss.ancom,detected_0.7 == TRUE)
sp.ss.ancom.sub<-subset(sp.ss.ancom,detected_0.7 == TRUE)
sp.sn.ancom.sub<-subset(sp.sn.ancom,detected_0.7 == TRUE)

### Combine DESeq and ANCOM II (0.7) tables-------------------------------------
sn.ss.combo<-inner_join(sn.ss.ancom.sub,sn.ss.flt,by="ASV")
sp.ss.combo<-inner_join(sp.ss.ancom.sub,sp.ss.flt,by="ASV")
sp.sn.combo<-inner_join(sp.sn.ancom.sub,sp.sn.flt,by="ASV")
