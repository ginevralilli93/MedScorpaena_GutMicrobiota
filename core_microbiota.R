# Core microbiota of Mediterranean Scorpionfishes as reported in the article Lilli et al., .
#
# The phyloseq object called in this script has been filtred previously for taxa_sums=0, mithocondria, cloroplasts and Archea, samples_sums<10000
# Comment: The taxaused are classified with the RDP classifier using the silva v138 database
# Species dataset (only BO) used for testing Species effect. Analyzed at the ASV level. 
# Extended dataset used for all the Geo analyses for each species separately --> datset including data from 3 studies = this study, Kormas et al., 2022, Ruiz-Rodriguez et al., 2020
# Extended dataset is analyzed only at the Genus level. 

# Author: Ginevra Lilli 

#Loading the packages that we will use
pkgs <- c("tidyverse", "phyloseq", "ggpubr", "ggplot2", 
          "vegan", "reshape2", "DESeq2","microbiome",
          "colorspace","picante","RColorBrewer","zCompositions","eulerr",
          "microbiomeutilities")
lapply(pkgs, require, character.only = TRUE)

# Generate the palettes----------------------------------------------------------
species_pal2<-c("royalblue4","orange2","yellow2")
mpa.palette.sn<-c("#FFCC00","grey60","grey2")
mpa.palette.sp.ss<-c("grey2","tan1")

# Functions----------------------------------------------------------------------
# remove the taxa with 0 counts in the new phy obj

nozero_taxa<-function(phy_rare){
  colsum.nozero<-(colSums(otu_table(phy_rare), na.rm=T) !=0)
  nozero1<-otu_table(phy_rare)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  phy_rare2<-phyloseq(otu_table(nozero1, taxa_are_rows = FALSE),
                      tax_table(tax_table(phy_rare)), taxa_names(taxa_names(phy_rare)),
                      sample_data(sample_data(phy_rare)),phy_tree(phy_tree(phy_rare)))
  
}

# Gather the data----------------------------------------------------------------

# Obtain the dada2+phyloseq object with additional samples from the 2 studies

# This is the Extended datset

sf<- get(load("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/merged_object/Rdata_obj/sf.filtered.object.with.tree.RData"))
meta.zone.location<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Data/metadata/125_scorpionfish_zone_by_loc.csv") 
meta.zone.location<-meta.zone.location[,-1]
rownames(meta.zone.location)<-meta.zone.location$FishID
sf<-phyloseq(otu_table(otu_table(sf), taxa_are_rows = FALSE),
             tax_table(tax_table(sf)), taxa_names(taxa_names(sf)),
             sample_data(meta.zone.location),phy_tree(phy_tree(sf)))

meta<-sample_data(sf) %>% data.frame()
meta %>% group_by(Species) %>% summarise(freq=n())


#1. Characterization of The Core Microbiota Of Three Species Of Scorpionfishes In The Mediterranean Sea"----------------------------------------------------------------------------

## 1.1. Species-dataset (BO)----------------------------------------------------

# Rarefy the data first to 20000 reads
sf.rare<-rarefy_even_depth(sf,rngseed = 123,sample.size = 20000, replace=FALSE)

# Subset the phyloseq object for the locations
species.dataset.BO<-subset_samples(sf.rare, MPA %in% 3) # 3 is the identifier of the BO region

species.dataset.BO<-nozero_taxa(species.dataset.BO)

# convert to relative abundances
View(otu_table(species.dataset.BO))
BO.gen<-tax_glom(species.dataset.BO,taxrank = "Genus", NArm = FALSE) 
BO.gen.ra<-transform_sample_counts(BO.gen,function(x){x / sum(x)})

#Make a list of the species
sf.sp <- unique(as.character(meta(BO.gen.ra)$Species))

#Write a for loop to go through each of the scorp species one by one and combine identified core taxa into a list
list_core.75.01.gen <- c() # an empty object to store information

for (n in sf.sp){ # for each variable n in Scorp.sp
  #print(paste0("Identifying Core Taxa for ", n))
  
  sf.rare.species.sub <- subset_samples(BO.gen.ra, Species == n) # Choose sample from Species by n
  
  core_m <- core_members(sf.rare.species.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001  
                         prevalence = 0.70)
  print(paste0("No. of core genera in ", n, " : ", length(core_m))) # print core taxa identified in each Species
  list_core.75.01.gen[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

# [1] "No. of core genera in Snotata : 20"
# [1] "No. of core genera in Sporcus : 17"
# [1] "No. of core genera in Sscrofa : 7"

p75.001.gen<-plot(venn(Alist_core.75.01.gen),
                  fills = species_pal2,
                  main = c("Genus level"),
                  quantities=TRUE,
                  cex.main=2)


list_core.sn.gen<-as.data.frame(list_core.75.01.gen$Snotata)
colnames(list_core.sn.gen)<-"ASV"



list_core.sp.gen<-as.data.frame(list_core.75.01.gen$Sporcus)
colnames(list_core.sp.gen)<-"ASV"

list_core.ss.gen<-as.data.frame(list_core.75.01.gen$Sscrofa)
colnames(list_core.ss.gen)<-"ASV"


# combine the taxonomy
tax.BO<-as.data.frame(tax_table(BO.gen))%>% rownames_to_column(var="ASV")
list_core.sn.gen<-inner_join(list_core.sn.gen,tax.BO,by="ASV")
list_core.sp.gen<-inner_join(list_core.sp.gen,tax.BO,by="ASV")
list_core.ss.gen<-inner_join(list_core.ss.gen,tax.BO,by="ASV")

# obtain the relative abundance of all the genera 

# subset for species 
sn.gen<-subset_samples(BO.gen.ra, Species == "Snotata")
sp.gen<-subset_samples(BO.gen.ra, Species == "Sporcus")
ss.gen<-subset_samples(BO.gen.ra, Species == "Sscrofa")

# extract the asv table 
sn.gen.asvtb<-as.data.frame(t(otu_table(sn.gen))) %>% rownames_to_column(var="ASV")
sp.gen.asvtb<-as.data.frame(t(otu_table(sp.gen))) %>% rownames_to_column(var="ASV")
ss.gen.asvtb<-as.data.frame(t(otu_table(ss.gen))) %>% rownames_to_column(var="ASV")

# join the asv table with the core taxa tables by ASV 
sn.core.gen<-inner_join(sn.gen.asvtb,list_core.sn.gen,by="ASV")
sp.core.gen<-inner_join(sp.gen.asvtb,list_core.sp.gen,by="ASV")
ss.core.gen<-inner_join(ss.gen.asvtb,list_core.ss.gen,by="ASV")

# write the csv files
write.csv2(sn.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/core_genera_with_abund_SN_narm.csv")
write.csv2(sp.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/core_genera_with_abund_SP_narm.csv")
write.csv2(ss.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/core_genera_with_abund_SS_narm.csv")

# Obtain the relative abundance of the core Genera on the full community 
# In excel: obtain the mean value of each Genus, calculate the "Not core" genera by 
# substracting the SUM of the core genera to 1 in each sample. Calcuate the average of 
# the "Not core" genera across all the samples. Reorder the MEAN values by ascending order.

abund.sn<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/sn_core_genera_mean_abund_narm.csv")
abund.sp<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/sp_core_genera_mean_abund_narm.csv")
abund.ss<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/ss_core_genera_mean_abund_narm.csv")


# Test for the significancy of the difference in the relative abundance of the core genera
# I will use a csv file I generated selecting only the sum of the relative abundances of each sample (from the files "core_genera_with_abundSN/SS/SP")
total<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/compare_Total_Core_Rel_Abund_narm.csv")

shapiro.test(total$Total) #not normal

#kruskal wallis 

kruskal.test(total$Total, total$Species) #  chi-squared = 9.79, df = 2, p-value = 0.007

# Dunn's test

FSA::dunnTest(total$Total, total$Species)

ggplot(total, aes(Species,Total))+
  geom_boxplot()+
  ylab("Relative abundance of core genera")+
  theme(axis.text.x = element_text(face="italic"))


# Find the shared and unique taxa 

shared.BO<-inner_join(abund.sn,abund.ss, by ="ASV") # 6 
shared.BO.sn.sp<-inner_join(abund.sn,abund.sp, by ="ASV")%>% subset(!ASV %in% shared.BO$ASV) %>% subset(!ASV %in% "Not core") # 8 
unique.sn<-subset(abund.sn, !ASV %in% c(shared.BO$ASV,shared.BO.sn.sp$ASV))%>% subset(!ASV %in% "Not core") # 6
unique.sp<-subset(abund.sp, !ASV %in% c(shared.BO$ASV,shared.BO.sn.sp$ASV))%>% subset(!ASV %in% "Not core") # 3
unique.ss<-subset(abund.ss, !ASV %in% c(shared.BO$ASV,shared.BO.sn.sp$ASV))
unique.ss<-subset(unique.ss,!ASV %in% "Not core " ) # 1

# write the csv files
write.csv2(shared.BO, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/shared_ALL.csv")
write.csv2(shared.BO.sn.sp, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/shared_SN_SP.csv")
write.csv2(unique.sn, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/unique_SN.csv")
write.csv2(unique.sp, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/unique_SP.csv")
write.csv2(unique.ss, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BO/unique_SS.csv")


## 1.2 Extended dataset for S.notata (BA,CR,CR)-----------------------------------

# Subset the phyloseq object for the locations
extended.dataset.SN<-subset_samples(sf.rare, MPA %in% c(1,2,3) & Species %in% "Snotata" & !FishID %in% "Microfish_3") # 3 is the identifier of the BO region

# Remove the txa with zero counts
extended.dataset.SN<-nozero_taxa(extended.dataset.SN)

# convert to relative abundances
EXT.gen<-tax_glom(extended.dataset.SN,taxrank = "Genus", NArm = FALSE) 
EXT.gen.ra<-transform_sample_counts(EXT.gen,function(x){x / sum(x)})

#Make a list of the species
sf.sn.EXT <- unique(as.character(meta(EXT.gen.ra)$MPA))

#Write a for loop to go through each of the scorp species one by one and combine identified core taxa into a list
list_core.75.01.genEXT <- c() # an empty object to store information

for (n in sf.sn.EXT){ # for each variable n in Scorp.sp
  #print(paste0("Identifying Core Taxa for ", n))
  
  sf.rare.geo.sub <- subset_samples(EXT.gen.ra, MPA == n) # Choose sample from Species by n
  
  core_m <- core_members(sf.rare.geo.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001  
                         prevalence = 0.75)
  print(paste0("No. of core genera in ", n, " : ", length(core_m))) # print core taxa identified in each Species
  list_core.75.01.genEXT[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

mpa.palette.sn<-c("#FFCC00","grey40","grey78")

p75.001.genEXT<-plot(venn(list_core.75.01.genEXT),
                  fills = mpa.palette.sn,
                  main = c("Genus level"),
                  quantities=TRUE,
                  cex.main=2)


list_core.BO.genEXT<-as.data.frame(list_core.75.01.genEXT$`3`)
colnames(list_core.BO.genEXT)<-"ASV"


list_core.CR.genEXT<-as.data.frame(list_core.75.01.genEXT$`2`)
colnames(list_core.CR.genEXT)<-"ASV"

list_core.BA.genEXT<-as.data.frame(list_core.75.01.genEXT$`1`)
colnames(list_core.BA.genEXT)<-"ASV"


# combine the taxonomy
tax.EXT<-as.data.frame(tax_table(EXT.gen))%>% rownames_to_column(var="ASV")
list_core.BO.genEXT<-inner_join(list_core.BO.genEXT,tax.EXT,by="ASV")
list_core.CR.genEXT<-inner_join(list_core.CR.genEXT,tax.EXT,by="ASV")
list_core.BA.genEXT<-inner_join(list_core.BA.genEXT,tax.EXT,by="ASV")

# obtain the relative abundance of all the genera 

# subset for Geo location  
BO.gen<-subset_samples(EXT.gen.ra, MPA == 3)
CR.gen<-subset_samples(EXT.gen.ra, MPA == 2)
BA.gen<-subset_samples(EXT.gen.ra, MPA == 1)

# extract the asv table 
BO.gen.asvtb<-as.data.frame(t(otu_table(BO.gen))) %>% rownames_to_column(var="ASV")
CR.gen.asvtb<-as.data.frame(t(otu_table(CR.gen))) %>% rownames_to_column(var="ASV")
BA.gen.asvtb<-as.data.frame(t(otu_table(BA.gen))) %>% rownames_to_column(var="ASV")

# join the asv table with the core taxa tables by ASV 
BO.core.gen<-inner_join(BO.gen.asvtb,list_core.BO.genEXT,by="ASV")
CR.core.gen<-inner_join(CR.gen.asvtb,list_core.CR.genEXT,by="ASV")
BA.core.gen<-inner_join(BA.gen.asvtb,list_core.BA.genEXT,by="ASV")

# write the csv files
write.csv2(BO.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/core_genera_with_abund_BO_narm.csv")
write.csv2(CR.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/core_genera_with_abund_CR_narm.csv")
write.csv2(BA.core.gen, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/core_genera_with_abund_BA_narm.csv")

# Obtain the relative abundance of the core Genera on the full community 
# In excel: obtain the mean value of each Genus, calculate the "Not core" genera by 
# substracting the SUM of the core genera to 1 in each sample. Calcuate the average of 
# the "Not core" genera across all the samples. Reorder the MEAN values by ascending order.

abund.BO<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/BO_core_genera_mean_abund_narm.csv")
abund.CR<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/CR_core_genera_mean_abund_narm.csv")
abund.BA<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/BA_core_genera_mean_abund_narm.csv")

# Find the shared and unique taxa 

shared.MED.SN<-inner_join(abund.BO,abund.BA, by ="ASV") # 7
shared.BO.CR<-inner_join(abund.BO,abund.CR, by ="ASV")%>% subset(!ASV %in% shared.MED.SN$ASV) %>% subset(!ASV %in% "Not core ") # 8 
shared.BA.CR<-inner_join(abund.BA,abund.CR, by ="ASV")%>% subset(!ASV %in% shared.MED.SN$ASV) %>% subset(!ASV %in% "Not core ") # 8 

unique.BO<-subset(abund.BO, !ASV %in% c(shared.MED.SN$ASV,shared.BO.CR$ASV))%>% subset(!ASV %in% "Not core ") # 9
unique.CR<-subset(abund.CR, !ASV %in% c(shared.MED.SN$ASV,shared.BA.CR$ASV,shared.BO.CR$ASV))%>% subset(!ASV %in% "Not core ") # 

# write the csv files
write.csv2(shared.MED.SN, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/shared_MED.csv")
write.csv2(shared.BO.CR, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/shared_BO_CR.csv")
write.csv2(shared.BA.CR, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/shared_BA_CR.csv")
write.csv2(unique.BO, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/unique_BO.csv")
write.csv2(unique.CR, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/BA_CR_BO_MPA_comparison/unique_CR.csv")


##1.3. Extended-dataset for S. porcus and S. scrofa (GY)---------------------------------------------------------------------------- 
# Rarefy the data first to 18025 reads to avoid to loose the samples from GY
sf.rareGy<-rarefy_even_depth(sf,rngseed = 123,sample.size = 18025, replace=FALSE)

# Subset the phyloseq object for the locations
species.dataset.GY<-subset_samples(sf.rareGy, MPA %in% 4) # 4 is the identifier of the GY region

# Remove the taxa with zero counts
species.dataset.GY<-nozero_taxa(species.dataset.GY)

# convert to relative abundances
GY.gen<-tax_glom(species.dataset.GY,taxrank = "Genus", NArm = FALSE) 
GY.gen.ra<-transform_sample_counts(GY.gen,function(x){x / sum(x)})

#Make a list of the species
sf.sp.Gy <- unique(as.character(meta(GY.gen.ra)$Species))

#Write a for loop to go through each of the scorp species one by one and combine identified core taxa into a list
Alist_core.75.01.genGY <- c() # an empty object to store information

for (n in sf.sp.Gy){ # for each variable n in Scorp.sp
  #print(paste0("Identifying Core Taxa for ", n))
  
  sf.rare.species.sub <- subset_samples(GY.gen.ra, Species == n) # Choose sample from Species by n
  
  core_m <- core_members(sf.rare.species.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001  
                         prevalence = 0.70)
  print(paste0("No. of core genera in ", n, " : ", length(core_m))) # print core taxa identified in each Species
  Alist_core.75.01.genGY[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

# [1] "No. of core genera in Sporcus : 14"
# [1] "No. of core genera in Sscrofa : 13"

Ap75.001.genGY<-plot(venn(Alist_core.75.01.genGY),
                  fills = species_pal2,
                  main = c("Genus level"),
                  quantities=TRUE,
                  cex.main=2)


list_core.sp.genGy<-as.data.frame(list_core.75.01.genGY$Sporcus)
colnames(list_core.sp.genGy)<-"ASV"

list_core.ss.genGy<-as.data.frame(list_core.75.01.genGY$Sscrofa)
colnames(list_core.ss.genGy)<-"ASV"


# combine the taxonomy
tax.GY<-as.data.frame(tax_table(GY.gen))%>% rownames_to_column(var="ASV")
list_core.sp.genGy<-inner_join(list_core.sp.genGy,tax.GY,by="ASV")
list_core.ss.genGy<-inner_join(list_core.ss.genGy,tax.GY,by="ASV")

# obtain the relative abundance of all the genera 

# subset for species 
sp.genGy<-subset_samples(GY.gen.ra, Species == "Sporcus")
ss.genGy<-subset_samples(GY.gen.ra, Species == "Sscrofa")

# extract the asv table 
sp.gen.asvtbGy<-as.data.frame(t(otu_table(sp.genGy))) %>% rownames_to_column(var="ASV")
ss.gen.asvtbGy<-as.data.frame(t(otu_table(ss.genGy))) %>% rownames_to_column(var="ASV")

# join the asv table with the core taxa tables by ASV 
sp.core.genGy<-inner_join(sp.gen.asvtbGy,list_core.sp.genGy,by="ASV")
ss.core.genGy<-inner_join(ss.gen.asvtbGy,list_core.ss.genGy,by="ASV")

# write the csv files
write.csv2(sp.core.genGy, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/core_genera_with_abund_SP_narm.csv")
write.csv2(ss.core.genGy, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/core_genera_with_abund_SS_narm.csv")

# Obtain the relative abundance of the core Genera on the full community 
# In excel: obtain the mean value of each Genus, calculate the "Not core" genera by 
# substracting the SUM of the core genera to 1 in each sample. Calcuate the average of 
# the "Not core" genera across all the samples. Reorder the MEAN values by ascending order.

abund.spGy<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/sp_core_genera_mean_abund_narm.csv")
abund.ssGy<-read.csv2("C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/ss_core_genera_mean_abund_narm.csv")


# Find the shared and unique taxa 

shared.Gy<-inner_join(abund.spGy,abund.ssGy, by ="ASV") # 6 
unique.spGy<-subset(abund.spGy, !ASV %in% c(shared.Gy$ASV))%>% subset(!ASV %in% "Not Core") # 3
unique.ssGy<-subset(abund.ssGy, !ASV %in% c(shared.Gy$ASV)) %>% subset(!ASV %in%  "Not core " ) # 1

# write the csv files
write.csv2(shared.Gy, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/shared_ALL.csv")
write.csv2(unique.spGy, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/unique_SP.csv")
write.csv2(unique.ssGy, "C:/Users/ginev/Dropbox/METRODIVER/WP3 MicroBiology/Analyses/10012022_16S/Output/outputs_comparison_external_data/core_microbiota_Second_Check/GY/unique_SS.csv")
















