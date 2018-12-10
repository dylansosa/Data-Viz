#####################################
# @author Dylan Sosa
# Dr. Rigdon
# Data Vizualiztion FA18
# Vizualizing Molecular Phylogenetics 
#####################################

# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# devtools::install_github("GuangchuangYu/ggtree")
# devtools::install_github('LaCroixColoR','johannesbjork')
require(ggplot2)
require(LaCroixColoR)
require(ggtree)

# load data 
setwd('/Users/dylansosa/Documents/SLU/5.3/Data Visualization/Programs_Data/project')
globins <- read.tree('vertebrateGlobins.tree')

# data exploration, partitioning 
globins$tip.label
groupInfo <- split(globins$tip.label, gsub("_[A-Z][a-z]+_.+", "",globins$tip.label))
groupInfo # the different globins 
globinTypes <- groupOTU(globins, groupInfo) # grouped based on globin type 

# re-root at outgroup
r <- root(globinTypes, node = 218, edgelabel = TRUE)

# draw phylogenetic tree
# 9 types of proteins from the Globin superfamily 
# here we can visualize the evolutionary history of these proteins in vertebrate species 
ggtree(r, aes(color=group)) + 
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
  geom_tiplab(size=1.4) + 
  geom_cladelabel(node=201, label="Beta Hemoglobin", offset = 0.7) + 
  geom_cladelabel(node=186, label="Alpha Hemoglobin", offset = 0.7) +
  geom_cladelabel(node=154, label="Cytoglobin", offset = 0.7) +
  geom_cladelabel(node=163, label="Agnathan Hemoglobin", offset = 0.7) +
  geom_cladelabel(node=174, label="Globin-Y", offset = 0.7) +
  geom_cladelabel(node=171, label="Globin-E", offset = 0.7) +
  geom_cladelabel(node=224, label="Neuroglobin", offset = 0.7) +
  geom_cladelabel(node=256, label="AndroGb", offset = 0.05) +
  geom_cladelabel(node=219, label="Outgroup", offset = 0.7, extend = 1.8) +
  geom_cladelabel(node=141, label="Myoglobin", offset = 0.7, extend = 2.1) +
  scale_color_manual(labels = groupInfo, 
                     values = lacroix_palette("Pamplemousse",type = "continuous", n=length(groupInfo)+1))

# load teneurin data 
teneurins <- read.tree('teneurins.tree')

# data exploration, partitioning 
teneurins$tip.label
teneuringroupInfo <- split(teneurins$tip.label, gsub("[A-Z][A-Z].+", "",teneurins$tip.label))
teneuringroupInfo # different species
groups <- groupOTU(teneurins, teneuringroupInfo) # species groups

# re-root with bacterial clade
rt <- root(teneurins, node = 152, edgelabel = TRUE)

# visualizing horizontal gene transfer of 
# polymorphic toxin proteins to eukaryotes from bacteria
# four copies of teneruins in vertebrates 
ggtree(rt)+
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  geom_tiplab(size=2) +
  geom_cladelabel(node=93, label="Teneurin 1", offset = 0.4) +
  geom_hilight(93,fill ='green') +
  geom_cladelabel(node=89, label="Teneurin 4", offset = 0.4) +
  geom_hilight(89,fill ='orange') +
  geom_cladelabel(node=98, label="Teneurin 2", offset = 0.4) +
  geom_hilight(98,fill ='blue') +
  geom_cladelabel(node=103, label="Teneurin 3", offset = 0.4) +
  geom_hilight(103,fill ='yellow') +
  geom_cladelabel(node=110, label="Basal Euk", offset = 0.4) +
  geom_hilight(110,fill ='brown') +
  geom_cladelabel(node=85, label="Coelomates") +
  geom_hilight(85,fill ='purple') 


# radial layout view to visualize differently
ggtree(teneurins,layout = 'daylight') +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  geom_tiplab(size=1.5) + 
  geom_hilight_encircle(node = 88) +
  geom_hilight_encircle(node = 92) +
  geom_hilight_encircle(node = 97) +
  geom_hilight_encircle(node = 102) +
  geom_hilight_encircle(node = 111, fill = 'darkgreen') +
  geom_hilight_encircle(node = 109, fill = 'brown') 

