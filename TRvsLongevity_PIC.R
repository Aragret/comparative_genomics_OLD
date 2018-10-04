rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

# user = 'Kostja'
user = 'Alina'

library(ape)

if (user == 'Alina'){
  setwd('/home/aragret/Alina/COMPARATIVE_GENOMICS/')
  CHOR <- as.data.frame(read.table("Derived_data/MitGenomics.txt", header = TRUE, sep = '\t'))
  tree <- read.tree("mtalign.aln.treefile.rooted")
}

library(ade4)

data = CHOR[which(as.character(CHOR$Species) %in% tree$tip.label),]

library(nlme)

df_vec <- as.character(CHOR$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

row.names(data) = data$Species

tree2 <- drop.tip(tree, b)


data$TRCoverage = data$REP.LengthOfTandemRepeats / data$GenomeLength
taxon_vec = unique(data$TAXON)

features <- c('TRCoverage', 'ECO.Maximum.longevity..yrs.')

row.names(data) = data$Species

one_line = c()
for(k in taxon_vec){
  TempData <- data[data$TAXON == k, features]
  for (i in 1:length(colnames(TempData))){
    TempData = TempData[!is.na(TempData[,i]),]
  }
  temp_diff = setdiff(tree2$tip.label, row.names(TempData))
  temp_tree <- drop.tip(tree2, temp_diff)
  
  contrasts <- apply(TempData, 2, pic, temp_tree)
  
  a = cor.test(log(contrasts[,1], base=exp(1)), log(contrasts[,2], base=exp(1)), method='spearman')
  one_line = rbind(one_line, c(features, k, nrow(TempData), a$estimate, a$p.value))
}

contrast_table <- as.data.frame(one_line)