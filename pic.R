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

fit2 <- gls(scale(GenomeLength) ~ 0 + scale(GCCont), data=data, correlation=corBrownian(phy=tree2))

str(fit2)
summary(fit2)

cs <- as.data.frame(summary(fit2)$tTable)
cs$p
# res <- resid(fit2, type="n")

##########################################################

names(CHOR)

features <- c('REP.NumberOfTandemRepeats', 'REP.LengthOfTandemRepeatsWithoutOverlaps', "REP.SymmRepNumber", "REP.InvRepNumber", "REP.DirRepNumber",
              "REP.ComplRepNumber", "GCCont", "GenomeLength", "GenomeWideGCSkew", "ECO.Female.maturity..days.",
              "ECO.Metabolic.rate..W.", "ECO.Body.mass..g.", 'ECO.Maximum.longevity..yrs.')

taxon_vec <- as.character(unique(CHOR$TAXON))

col_numbers = which(names(data) %in% features)

i = 6
j = 60
k = 1

temp_data = data[data$TAXON == taxon[k], c(i,j)]
temp_data = temp_data[(!is.na(temp_data[,1])),]
temp_data = temp_data[(!is.na(temp_data[,2])),]

temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
temp_tree <- drop.tip(tree, temp_diff)
colnames(temp_data) = c('features_1', 'features_2')

fit2 <- gls(features_1 ~ features_2, data=temp_data, correlation=corBrownian(phy=temp_tree))

### save species number! 

#####################################################################

VecOfColNumbers = which(names(data) %in% features)

data = CHOR[which(as.character(CHOR$Species) %in% tree$tip.label),]
one_line <- c()

for (i in VecOfColNumbers) 
{
  for (j in VecOfColNumbers)
  {
    if (j > i)
    {
      
      for (k in 0:length(taxon))
        {
        if (k == 0){
          TempData = data[,c(i,j)]
        }
        else {
          TempData = data[data$TAXON == taxon[k], c(i,j)]
        }
        TempData = TempData[!is.na(TempData[1]),]; TempData = TempData[!is.na(TempData[2]),]; 
        names(TempData) = c('FirstVar','SecondVar') 
        if (nrow(TempData) > 10){
          TempDiff <- setdiff(tree2$tip.label, row.names(TempData)); print(length(TempDiff))
          TempTree <- drop.tip(tree2, TempDiff); print(length(TempTree$tip.label)) # 191
          TempRes <- gls(FirstVar ~ SecondVar,  data = TempData, correlation=corBrownian(phy=TempTree))
          summary(TempRes)
          # TempRes$coefficients
          
          cs <- as.data.frame(summary(TempRes)$tTable)
          cs$p
          if(k == 0){
            one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), 'All', TempRes$coefficients, cs$p))
            next
          }
          one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), taxon[k], TempRes$coefficients, cs$p))
        }
      }
    }
  }
}

Result <- as.data.frame(one_line)
names(Result) <- c('FirstVar', 'SecondVar', 'SpeciesNumber', 'TAXON', 'Intersept',
                   'Slope', 'InterseptPValue', 'SlopePValue')

write.table(Result, file = 'Derived_data/gls.txt', sep = '\t', quote = FALSE, row.names = FALSE)

######################################################################

temp_data = data[data$TAXON == taxon[k], c(i,j)]
temp_data = temp_data[(!is.na(temp_data[,1])),]
temp_data = temp_data[(!is.na(temp_data[,2])),]

temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
temp_tree <- drop.tip(tree, temp_diff)
colnames(temp_data) = c('features_1', 'features_2')


features_1 = 'GCCont'
features_2 = 'REP.NumberOfTandemRepeats'
temp_data = data[, c(features_1, features_2)]
colnames(temp_data) = c('features_1', 'features_2')
temp_data = temp_data[!is.na(temp_data$features_2), ]
temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
temp_tree <- drop.tip(tree2, temp_diff)
pglsModelLambda <- gls(features_1 ~ features_2, correlation = corPagel(1, phy = temp_tree,
                                                                       fixed = FALSE), data = temp_data, method = "ML")
features_1 = 'GCCont'
features_2 = 'ECO.Maximum.longevity..yrs'
temp_data = data[, c(features_1, features_2)]
colnames(temp_data) = c('features_1', 'features_2')
temp_data = temp_data[!is.na(temp_data$features_2), ]
temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
temp_tree <- drop.tip(tree2, temp_diff)
pglsModelLambda <- gls(features_1 ~ features_2, correlation = corMartins(1, phy = temp_tree,
                                                                       fixed = FALSE), data = temp_data, method = "ML")




cs <- as.data.frame(summary(fit)$tTable)
cs$p

#########################################################################

VecOfColNumbers = which(names(data) %in% features)

data = CHOR[which(as.character(CHOR$Species) %in% tree$tip.label),]
one_line <- c()
taxon <- as.character(unique(CHOR$TAXON))
row.names(data) = data$Species

for (i in VecOfColNumbers) 
{
  for (j in VecOfColNumbers)
  {
    if (j > i)
    {
      
      for (k in 0:length(taxon))
      {
        if (k == 0){
          TempData = data[,c(i,j)]
        }
        else {
          TempData = data[data$TAXON == taxon[k], c(i,j)]
        }
        TempData = TempData[!is.na(TempData[1]),]; TempData = TempData[!is.na(TempData[2]),]; 
        names(TempData) = c('FirstVar','SecondVar') 
        if (nrow(TempData) > 10){
          TempDiff <- setdiff(tree2$tip.label, row.names(TempData)); print(length(TempDiff))
          TempTree <- drop.tip(tree2, TempDiff); print(length(TempTree$tip.label)) # 191
          TempRes <- gls(scale(FirstVar) ~ scale(SecondVar),  data = TempData, correlation=corBrownian(phy=TempTree))
          summary(TempRes)
          # TempRes$coefficients
          
          cs <- as.data.frame(summary(TempRes)$tTable)
          cs$p
          if(k == 0){
            one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), 'All', TempRes$coefficients, cs$p))
            next
          }
          one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), taxon[k], TempRes$coefficients, cs$p))
        }
      }
    }
  }
}

Scaled_result <- as.data.frame(one_line)
names(Scaled_result) <- c('FirstVar', 'SecondVar', 'SpeciesNumber', 'TAXON', 'Intersept',
                   'Slope', 'InterseptPValue', 'SlopePValue')

#############################################################
### Pagel

VecOfColNumbers = which(names(data) %in% features)
one_line <- c(NULL)
count = 0

for (i in VecOfColNumbers) 
{
  for (j in VecOfColNumbers)
  {
    if (j > i)
    {
      
      for (k in 0:length(taxon))
      {
        # i = 21; j = 52; k = 4
        if (k == 0){
          # TempData = data[,c(i,j)]
          next
        }
        else {
          TempData = data[data$TAXON == taxon[k], c(i,j)]
        }
        TempData = TempData[!is.na(TempData[1]),]; TempData = TempData[!is.na(TempData[2]),]; 
        names(TempData) = c('FirstVar','SecondVar') 
        if (nrow(TempData) > 10){
          TempDiff <- setdiff(tree2$tip.label, row.names(TempData)); print(length(TempDiff))
          TempTree <- drop.tip(tree2, TempDiff); print(length(TempTree$tip.label)) # 191
          TempRes <- try(gls(scale(FirstVar) ~ scale(SecondVar),  data = TempData, correlation=corPagel(0.8, phy=TempTree), method = 'ML'), silent = TRUE)
          a = summary(TempRes)
          # TempRes$coefficients
          if (a[2] != 'try-error'){
            cs <- as.data.frame(summary(TempRes)$tTable)
            count = count + 1
            if(k == 0){
              one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), 'All', TempRes$coefficients, cs$p,
                                            a$modelStruct, a$logLik))
              next
            }
            one_line <- rbind(one_line, c(colnames(data[,c(i,j)]), nrow(TempData), taxon[k], TempRes$coefficients, cs$p,
                                          a$modelStruct, a$logLik))
          }
        }
      }
    }
  }
}

first_result = as.data.frame(one_line)
names(first_result) = c('Feature_1', 'Feature_2', 'SpeciesNumber', 'TAXON',
                        'Intersept', 'Slope', 'InterseptPValue', 'SlopePValue',
                        'Lambda', 'LogLik')

# write.table(first_result, file = 'PagelGLS.csv', sep = '\t', quote = FALSE, row.names = FALSE) # 
library(data.table)
fwrite(first_result, file ="PagelGLSScaled.csv")

features_1 = 'REP.NumberOfTandemRepeats'
features_2 = 'ECO.Female.maturity..days.'
temp_data = data[data$TAXON == 'Reptilia', c(features_1, features_2)]
colnames(temp_data) = c('features_1', 'features_2')
temp_data = temp_data[!is.na(temp_data$features_2), ]
temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
temp_tree <- drop.tip(tree2, temp_diff)
pglsModelLambda <- gls(features_1 ~ features_2, correlation = corMartins(1, phy = temp_tree,
                                                                         fixed = FALSE), data = temp_data, method = "ML")
try(gls(features_1 ~ features_2, correlation = corMartins(1, phy = temp_tree,
                                                          fixed = FALSE), data = temp_data, method = "ML"),
    silent = TRUE)

#########################################################################
### multiple models

features <- c('REP.NumberOfTandemRepeats', 'REP.DirRepLength', 'REP.SymmRepLength',
              'REP.ComplRepLength', 'REP.InvRepLength', 'GCCont', 'REP.LengthOfTandemRepeatsWithoutOverlaps')

row.names(data) = data$Species
TempData <- data[data$TAXON == 'Mammalia', c('GenomeLength', features)]
for (i in 1:length(colnames(TempData))){
  TempData = TempData[!is.na(TempData[,i]),]
}
temp_diff = setdiff(tree2$tip.label, row.names(TempData))
temp_tree <- drop.tip(tree2, temp_diff)

fit = gls(scale(GenomeLength) ~ 0 + (scale(REP.DirRepNumber) +
            scale(REP.ComplRepNumber)) * scale(GCCont) +
            scale(REP.NumberOfTandemRepeats) + scale(REP.SymmRepNumber),
          data = TempData, correlation = corPagel(0.3, phy = temp_tree), method = 'ML')

fit = gls(scale(GenomeLength) ~ 0 + scale(REP.DirRepNumber) + scale(REP.ComplRepNumber) +
            scale(REP.NumberOfTandemRepeats) + scale(REP.SymmRepNumber) + scale(REP.InvRepNumber) +
            scale(GCCont),
          data = TempData, correlation = corPagel(0.3, phy = temp_tree), method = 'ML')



lambda <- seq(0, 1, length.out = 10)
for (i in lambda){
  fit = try(gls(GenomeLength ~ REP.NumberOfTandemRepeats + 
                  REP.DirRepNumber + REP.SymmRepNumber +
                  REP.ComplRepNumber + REP.InvRepNumber,
                data = TempData, correlation = corPagel(i, phy = temp_tree), method = 'ML'))
  a = summary(fit)
  if (a[2] != 'try-error'){
    print(c(a$modelStruct, a$logLik))
  }
  next
}


fit = try(gls(log(GenomeLength, base = exp(1)) ~ log(REP.NumberOfTandemRepeats, base = exp(1)) + 
            log(REP.DirRepNumber, base = exp(1)) + log(REP.SymmRepNumber, base = exp(1)) +
            log(REP.ComplRepNumber, base = exp(1)) + log(REP.InvRepNumber, base = exp(1)),
          data = TempData, correlation = corPagel(0.3, phy = temp_tree), method = 'ML'))


###################################################################################
### residuals

fit = gls(scale(GenomeLength) ~ 0 + scale(REP.DirRepNumber) +
            scale(REP.ComplRepNumber) + scale(REP.NumberOfTandemRepeats)
          + scale(REP.SymmRepNumber),
          data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')

res <- residuals(fit)

summary(res)

a = gls(res ~ 0 + GCCont, data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')
summary(a)

cor.test(res, TempData$GCCont, method = 'spearman')

a = gls(res ~ 0 + ECO.Female.maturity..days., data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')
summary(a)

##########################################################
### repeat length

fit = gls(scale(GenomeLength) ~ 0 + scale(REP.DirRepLength) +
            scale(REP.ComplRepLength) + scale(REP.LengthOfTandemRepeatsWithoutOverlaps)
          + scale(REP.SymmRepLength),
          data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')

res <- resid(fit, type='n')

a = gls(res ~ 0 + GCCont, data = TempData, correlation = corPagel(0.3, phy = temp_tree), method = 'ML')
summary(a)

cor.test(res, TempData$GCCont, method = 'spearman')

a = gls(res ~ 0 + ECO.Female.maturity..days., data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')
summary(a)

################################################################################################################################

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'REP.DirRepLength', 'REP.SymmRepLength', 'REP.ComplRepLength',
              'REP.InvRepLength')

row.names(data) = data$Species

one_line <- c()
for (tax in taxon_vec){
  TempData <- data[data$TAXON == tax, features]
  for (i in 1:length(colnames(TempData))){
    TempData = TempData[!is.na(TempData[,i]),]
  }
  temp_diff = setdiff(tree2$tip.label, row.names(TempData))
  temp_tree <- drop.tip(tree2, temp_diff)
  
  fit = gls(scale(GenomeLength) ~ 0 + scale(REP.DirRepLength) +
              scale(REP.ComplRepLength) + scale(REP.LengthOfTandemRepeatsWithoutOverlaps)
            + scale(REP.SymmRepLength) + scale(REP.InvRepLength),
            data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')
  a = summary(fit)
  cs <- as.data.frame(summary(fit)$tTable)
  one_line <- rbind(one_line, c(tax, a$coefficients, cs$p, a$modelStruct, a$logLik))
  
}

tab = as.data.frame(one_line)
names(tab) = c('TAXON', 'SlopeDirLength', 'SlopeComplLength', 'SlopeTandemLength',
               'SlopeSymmLength', 'SlopeInvLength', 'PValueDirLength', 'PValueComplLength', 
               'PValueTandemLength', 'PValueSymmLength', 'PValueInvLength',
               'Lambda', 'Loglik')

library(data.table)
fwrite(tab, file = '../Alina/COMPARATIVE_GENOMICS/Derived_data/GLvsRepeatsPagel.csv')

###########################################################################
### PIC

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'REP.DirRepLength', 'REP.SymmRepLength', 'REP.ComplRepLength',
              'REP.InvRepLength')

row.names(data) = data$Species

ContrastwingL <- pic(wingL, geotree)
ContrasttarsusL <- pic(tarsusL, geotree)

one_line <- c()
for (i in 1:length(features)){
  a = c(paste('PIC_', i, sep=''), pic(data[,features[i]], tree2))
  one_line = rbind(one_line, a)
}

TempData <- one_line[c(1, 2),-1]
first_vec = as.numeric(one_line[1,-1])

lm(as.vector(one_line[1,-1]) ~ as.vector(one_line[2,-1]))

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'REP.DirRepLength', 'REP.SymmRepLength', 'REP.ComplRepLength',
              'REP.InvRepLength')

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

  a = summary(lm(contrasts[,1] ~ 0 + contrasts[,2] +  contrasts[,3] +
                   contrasts[,4] + contrasts[,5] + contrasts[,6]))
  one_line = rbind(one_line, c(k, a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SlopeTandemLength', 'SlopeDirLength', 'SlopeSymmLength',
               'SlopeComplLength', 'SlopeInvLength', 'PValueTandLength', 'PValueDirLength', 
               'PValueSymmLength', 'PValueComplLength', 'PValueInvLength')

write.table(contrast_table, file = '../Alina/COMPARATIVE_GENOMICS/Derived_data/GLvsRepeatsPIC.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

#########################################################################
### TR and GC vs GL

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont')
taxon_vec <- c('Actinopterygii', 'Reptilia', 'Aves', 'Mammalia', 'Amphibia')

row.names(data) = data$Species

one_line <- c()
count = 0
for (tax in taxon_vec){
  TempData <- data[data$TAXON == tax, features]
  for (i in 1:length(colnames(TempData))){
    TempData = TempData[!is.na(TempData[,i]),]
  }
  temp_diff = setdiff(tree2$tip.label, row.names(TempData))
  temp_tree <- drop.tip(tree2, temp_diff)
  
  fit = try(gls(scale(GenomeLength) ~ 0 + scale(REP.LengthOfTandemRepeatsWithoutOverlaps)
            + scale(GCCont),
            data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML'))
  a = summary(fit)
  if(a[2] == 'try-error'){
    fit = gls(scale(GenomeLength) ~ 0 + scale(REP.LengthOfTandemRepeatsWithoutOverlaps)
                  + scale(GCCont),
                  data = TempData, correlation = corPagel(0.8, phy = temp_tree), method = 'ML')
  }
  cs <- as.data.frame(summary(fit)$tTable)
  one_line <- rbind(one_line, c(tax, a$coefficients, cs$p, a$modelStruct, a$logLik))
}

tab = as.data.frame(one_line)
names(tab) = c('TAXON', 'SlopeDirLength', 'SlopeComplLength', 'SlopeTandemLength',
               'SlopeSymmLength', 'SlopeInvLength', 'PValueDirLength', 'PValueComplLength', 
               'PValueTandemLength', 'PValueSymmLength', 'PValueInvLength',
               'Lambda', 'Loglik')

library(data.table)
fwrite(tab, file = '../Alina/COMPARATIVE_GENOMICS/Derived_data/GLvsRepeatsPagel.csv')


### with PIC

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont')

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
  # a = summary(lm(contrasts[,1] ~ 0 + contrasts[,2] + contrasts[,3]))
  a = summary(lm(scale(contrasts[,1]) ~ 0 + scale(contrasts[,2]) +  
                   scale(contrasts[,3])))
  one_line = rbind(one_line, c(k, a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SlopeTandemLength', 'SlopeGCCont', 'PValueTandLength', 
                          'PValueGCCont')

write.table(contrast_table, 'MultipleModels/ScaledGLvsTRplusGC.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

### 

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont')

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
  
  a = summary(lm(scale(contrasts[,1]) ~ 0 + scale(contrasts[,2]) *  scale(contrasts[,3])))
  one_line = rbind(one_line, c(k, a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SlopeTandemLength', 'SlopeGCCont', 'SlopeInteraction', 
                          'PValueTandLength', 'PValueGCCont', 'PValueInteraction')

write.table(contrast_table, 'MultipleModels/ScaledGLvsTRandGCInteract.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

############################

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont', 'ECO.Female.maturity..days.')

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
  
  a = summary(lm(scale(contrasts[,1]) ~ 0 + scale(contrasts[,2]) + scale(contrasts[,3]) + 
                   scale(contrasts[,4])))
  one_line = rbind(one_line, c(k, nrow(TempData), a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SpeciesNumber', 'SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity')

write.table(contrast_table, 'MultipleModels/ScaledGLvsTRplusGCplusFM.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

###

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont', 'ECO.Female.maturity..days.')

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
  
  a = summary(lm(scale(contrasts[,1]) ~ 0 + scale(contrasts[,2]) + scale(contrasts[,3])
                 * scale(contrasts[,4])))
  one_line = rbind(one_line, c(k, nrow(TempData), a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SpeciesNumber','SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'SlopeFemaleMaturityAndGC',
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity', 'PValueFemaleMaturityAndGC')

write.table(contrast_table, 'MultipleModels/ScaledGLvsTRplusGC*FM.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

###

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont', 'ECO.Female.maturity..days.')

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
  
  a = summary(lm(scale(contrasts[,1]) ~ 0 + (scale(contrasts[,2]) + 
                                               scale(contrasts[,3])) * 
                   scale(contrasts[,4])))
  one_line = rbind(one_line, c(k, nrow(TempData), a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SpeciesNumber','SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'SlopeTandemLengthAndGC', 'SlopeFemaleMaturityAndGC',
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity',
                          'PValueTandemLengthAndGC', 'PValueFemaleMaturityAndGC')

write.table(contrast_table, 'MultipleModels/ScaledGLvs(TRplusFM)*GC.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'GCCont', 'ECO.Female.maturity..days.')

######################################################################
### All pairwise PIC

row.names(data) = data$Species

features <- c('REP.LengthOfTandemRepeatsWithoutOverlaps', 'REP.DirRepLength',
              'REP.SymmRepNumber', 'REP.ComplRepNumber', 'REP.InvRepNumber',
              'GenomeLength', 'GCCont', 'GenomeWideGCSkew', 'ECO.Female.maturity..days.',
              'ECO.Maximum.longevity..yrs.')

one_line = c()
for(k in taxon_vec){
  TempData <- data[data$TAXON == k, features]
  for (i in 1:length(colnames(TempData))){
    TempData = TempData[!is.na(TempData[,i]),]
  }
  for(i in 1:length(features)){
    for(j in 1:length(features)){
      if(j > i){
        TempData <- data[data$TAXON == k, c(features[i], features[j])]
        TempData = TempData[!is.na(TempData[,1]),]; TempData = TempData[!is.na(TempData[,2]),]
        temp_diff = setdiff(tree2$tip.label, row.names(TempData))
        temp_tree <- drop.tip(tree2, temp_diff)
        contrasts <- apply(TempData, 2, pic, temp_tree)
        a = cor.test(scale(contrasts[,1]), scale(contrasts[,2]), method = 'spearman')
        one_line = rbind(one_line, c(features[i], features[j], k, nrow(TempData),
                                     as.character(a$estimate),
                                     as.character(a$p.value)))
      }
    }
  }
  
}

results = as.data.frame(one_line)
names(results) = c('FirstFeature', 'SecondReature', 'TAXON', 'SpeciesNumber',
                   'Rho', 'PValue')

write.table(results, file = 'Derived_data/PairwisePIC.txt', sep='\t',
            quote=FALSE, row.names = FALSE)

####################

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
names(contrast_table) = c('TAXON', 'SpeciesNumber','SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'SlopeTandemLengthAndGC', 'SlopeFemaleMaturityAndGC',
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity',
                          'PValueTandemLengthAndGC', 'PValueFemaleMaturityAndGC')

write.table(contrast_table, 'MultipleModels/ScaledGLvs(TRplusFM)*GC.txt', sep='\t', quote = FALSE,
            row.names = FALSE)

plot(log(contrasts[,1], base=exp(1)), log(contrasts[,2], base=exp(1)))



features <- c('REP.LengthOfTandemRepeatsWithoutOverlaps', 'REP.DirRepLength',
              'REP.SymmRepNumber', 'REP.ComplRepNumber', 'REP.InvRepNumber',
              'GenomeLength', 'GCCont', 'GenomeWideGCSkew')
TempData = data[,features]
contrasts <- apply(TempData, 2, pic, tree2)

features <- c('ECO.Female.maturity..days.', 'ECO.Maximum.longevity..yrs.')
TempData = data[,features]
TempFeature = TempData[!is.na(TempData[,1]),]
temp_diff = setdiff(tree2$tip.label, row.names(TempFeature))
temp_tree <- drop.tip(tree2, temp_diff)
a = cbind(contrasts, pic(TempFeature[,1], temp_tree))

TempFeature = TempData[!is.na(TempData[,2]),]
temp_diff = setdiff(tree2$tip.label, row.names(TempFeature))
temp_tree <- drop.tip(tree2, temp_diff)
a = cbind(a, pic(TempFeature[,2], temp_tree))

cor.test(a[,9], a[,10])

### loop for calculating contrasts for each feature
features <- c('REP.LengthOfTandemRepeatsWithoutOverlaps', 'REP.DirRepLength',
              'REP.SymmRepNumber', 'REP.ComplRepNumber', 'REP.InvRepNumber',
              'GenomeLength', 'GCCont', 'GenomeWideGCSkew', 'ECO.Female.maturity..days.',
              'ECO.Maximum.longevity..yrs.')

########################################################################################
### table with contrasts

library(phytools)
features <- c('REP.LengthOfTandemRepeatsWithoutOverlaps', 'REP.DirRepLength',
              'REP.SymmRepNumber', 'REP.ComplRepNumber', 'REP.InvRepNumber',
              'GenomeLength', 'GCCont', 'GenomeWideGCSkew', 'ECO.Female.maturity..days.',
              'ECO.Maximum.longevity..yrs.')
row.names(data) = data$Species

one_line <- c()

for(i in features){
  # i = 'ECO.Maximum.longevity..yrs.'
  # i = 'REP.LengthOfTandemRepeatsWithoutOverlaps'
  temp_data = data[, c('Species', i)]
  temp_data = temp_data[!is.na(temp_data[, i]),]
  temp_diff = setdiff(tree2$tip.label, row.names(temp_data))
  temp_tree <- drop.tip(tree2, temp_diff)
  contrasts <- pic(temp_data[, i], temp_tree)
  print(length(contrasts))
  #print(nodelabels(round(contrasts.var [,1], 3), adj = c(0, -0.5), frame="n"))
  one_line = rbind(one_line, c(i, paste(contrasts, collapse = ' ')))
}

one_line = as.data.frame(one_line)
names(one_line) = c('Features', 'Contrasts')

# write.table(one_line, x = 'contrasts_table.txt', sep = '\t', quote = FALSE, row.names = FALSE)
library(data.table)
fwrite(one_line, file = 'Derived_data/contrasts_table.csv')

##################################################
### create table for multiple analysis

data$TRCoverage = data$REP.LengthOfTandemRepeats / data$GenomeLength

for(i in 1:nrow(data)){
  if(data$REP.NumberOfTandemRepeats[i] == 0){
    data$TRPresence[i] = 0
  }
  else{data$TRPresence[i] = 1}
}

features <- c('GenomeLength', 'GCCont', 'TRPresence',
              'TRCoverage',  'ECO.Female.maturity..days.')
row.names(data) = data$Species

one_line = c()

TempData <- data[data$TAXON == 'Mammalia', features]
for (i in 1:length(colnames(TempData))){
  
  TempData = TempData[!is.na(TempData[,i]),]
}
  
temp_diff = setdiff(tree2$tip.label, row.names(TempData))
temp_tree <- drop.tip(tree2, temp_diff)
  
contrasts <- apply(TempData, 2, pic, temp_tree)
  
a = cor.test(log(contrasts[,1], base=exp(1)), log(contrasts[,2], base=exp(1)), method='spearman')
one_line = rbind(one_line, c(features, k, nrow(TempData), a$estimate, a$p.value))

contrast_table <- as.data.frame(contrasts)
names(contrast_table) = c('TAXON', 'SpeciesNumber','SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'SlopeTandemLengthAndGC', 'SlopeFemaleMaturityAndGC',
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity',
                          'PValueTandemLengthAndGC', 'PValueFemaleMaturityAndGC')

write.table(contrast_table, 'Derived_data/TableForMultipleModels', sep='\t', quote = FALSE,
            row.names = FALSE)


# write.table(one_line, x = 'contrasts_table.txt', sep = '\t', quote = FALSE, row.names = FALSE)
library(data.table)
fwrite(one_line, file = 'Derived_data/contrasts_table.csv')


###

features <- c('TRCoverage',  'ECO.Female.maturity..days.')
row.names(data) = data$Species

one_line = c()

TempData <- data[data$TAXON == 'Mammalia', features]
for (i in 1:length(colnames(TempData))){
  
  TempData = TempData[!is.na(TempData[,i]),]
}

temp_diff = setdiff(tree2$tip.label, row.names(TempData))
temp_tree <- drop.tip(tree2, temp_diff)

contrasts <- apply(TempData, 2, pic, temp_tree)

a = cor.test(log(contrasts[,1], base=exp(1)), log(contrasts[,2], base=exp(1)), method='spearman')
one_line = rbind(one_line, c(features, nrow(TempData), a$estimate, a$p.value))

contrast_table <- as.data.frame(contrasts)
names(contrast_table) = c('TAXON', 'SpeciesNumber','SlopeTandemLength', 'SlopeGCCont', 'SlopeFemaleMaturity', 
                          'SlopeTandemLengthAndGC', 'SlopeFemaleMaturityAndGC',
                          'PValueTandLength', 'PValueGCCont', 'PValueFemaleMaturity',
                          'PValueTandemLengthAndGC', 'PValueFemaleMaturityAndGC')

write.table(contrast_table, 'Derived_data/TableForMultipleModels', sep='\t', quote = FALSE,
            row.names = FALSE)


# write.table(one_line, x = 'contrasts_table.txt', sep = '\t', quote = FALSE, row.names = FALSE)
library(data.table)
fwrite(one_line, file = 'Derived_data/contrasts_table.csv')


library(ggplot2)
ggplot(data = contrast_table, aes(log(TRCoverage, base = exp(1)), 
                                  log(ECO.Female.maturity..days., base = exp(1)))) +
  geom_point()
