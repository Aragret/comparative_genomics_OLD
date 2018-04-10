####################################################################################
############ Mammalian Tandem Repeats with high Percent Match are more GC rich
####################################################################################

rm(list=ls(all=TRUE))

##### A: split all TRs into young (PM > 95%) and old (PM < 95%) and analyze mammalian species with existed both categories

user = 'Kostya'
if (user == 'Kostya')
{
MM <- read.table("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/3_RESULTS/MitGenomics.txt", header = TRUE, sep = '\t')  
YoungTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/YoungTandemRepeats.txt', header = TRUE, sep = '\t')
OldTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/OldTandemRepeats.txt', header = TRUE, sep = '\t')
}

MM = MM[MM$TAXON == 'Mammalia',]
VecOfMammalianSpecies = as.character(MM[,1]); length(VecOfMammalianSpecies)
names(MM)[1] = c('Species')

YoungAgg = aggregate(YoungTr$ATCont, by = list(YoungTr$Species), FUN = mean)
names(YoungAgg) = c('Species','YoungTrMeanAt')

OldAgg = aggregate(OldTr$ATCont, by = list(OldTr$Species), FUN = mean)
names(OldAgg) = c('Species','OldTrMeanAt')

if (user == 'Kostya')
{pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/4_FIGURES/TandemRepeatsPmVsAt.pdf')}

par(mfrow=c(2,2))
Tr = merge(YoungAgg,OldAgg, by = 'Species')
Tr = Tr[Tr$Species %in% VecOfMammalianSpecies,] 
nrow(Tr) # 566
hist(Tr$YoungTrMeanAt, breaks = 50, xlab = 'AT content of TRs with high Percent Match'); 
hist(Tr$OldTrMeanAt, breaks = 50, xlab = 'AT content of TRs with low Percent Match'); 

Tr$YoungAtMinusOldAt = Tr$YoungTrMeanAt - Tr$OldTrMeanAt
summary(Tr$YoungAtMinusOldAt)
hist(Tr$YoungAtMinusOldAt, breaks = 100, xlab = 'AtInHighPm-AtInLowPm')
abline(v = median(Tr$YoungAtMinusOldAt), col = 'red', lwd = 2)
wilcox.test(Tr$YoungAtMinusOldAt, mu = 0)
wilcox.test(Tr$YoungTrMeanAt,Tr$OldTrMeanAt, paired = TRUE)
t.test(Tr$YoungTrMeanAt,Tr$OldTrMeanAt, paired = TRUE)
plot(Tr$YoungTrMeanAt,Tr$OldTrMeanAt)
abline(a = 0, b= 1, col = 'red')

##### B: Analyse correlation withing each genome with at least two TRs: the higher PM the lower AT 

if (user == 'Kostya')
{
  MM <- read.table("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/3_RESULTS/MitGenomics.txt", header = TRUE, sep = '\t')  
  YoungTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/YoungTandemRepeats.txt', header = TRUE, sep = '\t')
  OldTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/OldTandemRepeats.txt', header = TRUE, sep = '\t')
}

Tr = rbind(YoungTr,OldTr)
Tr$TrLength = Tr$End - Tr$Start
Tr = Tr[Tr$Species %in% VecOfMammalianSpecies,] 

Tr$Number = 1
AGG = aggregate(Tr$Number, by = list(Tr$Species), FUN = sum)
summary(AGG$x) # .000   1.000   2.000   2.132   3.000  11.000 
VecOfSpeciesWithManyTr = as.character(AGG[AGG$x >= 2,]$Group.1); length(VecOfSpeciesWithManyTr)
Tr = Tr[Tr$Species %in% VecOfSpeciesWithManyTr,]
for (i in 1:length(VecOfSpeciesWithManyTr))
{ # i = 1
Temp = Tr[Tr$Species %in% VecOfSpeciesWithManyTr[i],]  
PmVsAtPval = as.numeric(cor.test(Temp$PercentMatches,Temp$ATCont, method = 'spearman')[3])
PmVsAtRho = as.numeric(cor.test(Temp$PercentMatches,Temp$ATCont, method = 'spearman')[4])
OneLine = data.frame(VecOfSpeciesWithManyTr[i],PmVsAtPval,PmVsAtRho)
if (i == 1) {Final=OneLine}
if (i >  1) {Final=rbind(Final,OneLine)}
}

# high Pm -> low AT
summary(Final$PmVsAtRho) # negative
summary(Final[Final$PmVsAtPval < 0.1,]$PmVsAtRho) # negative
hist(Final[Final$PmVsAtPval < 0.1,]$PmVsAtRho, breaks = 50, xlab = 'Rho with p < 0.1')
plot(Final$PmVsAtRho, -log2(Final$PmVsAtPval), xlab = 'Rho', ylab = '-log2(p)')
abline(h = -log2(0.1), col = 'red', lwd = 3)

##### C: Between Species approach = all TRs with different PM analyse in one mix

if (user == 'Kostya')
{
  MM <- read.table("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/3_RESULTS/MitGenomics.txt", header = TRUE, sep = '\t')  
  YoungTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/YoungTandemRepeats.txt', header = TRUE, sep = '\t')
  OldTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/OldTandemRepeats.txt', header = TRUE, sep = '\t')
}

Tr = rbind(YoungTr,OldTr)
Tr$TrLength = Tr$End - Tr$Start
Tr = Tr[Tr$Species %in% VecOfMammalianSpecies,] 

summary(Tr$PercentMatches)
boxplot(Tr[Tr$PercentMatches < quantile(Tr$PercentMatches,0.5),]$ATCont,Tr[Tr$PercentMatches >= quantile(Tr$PercentMatches,0.5),]$ATCont, notch = TRUE, names = c('old','young'), ylab = 'AT Content', main = 'interspecies comparison')


##### D: more accurate between Species approach = compare averages

if (user == 'Kostya')
{
  MM <- read.table("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/3_RESULTS/MitGenomics.txt", header = TRUE, sep = '\t')  
  YoungTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/YoungTandemRepeats.txt', header = TRUE, sep = '\t')
  OldTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/OldTandemRepeats.txt', header = TRUE, sep = '\t')
}

Tr = rbind(YoungTr,OldTr)
Tr$TrLength = Tr$End - Tr$Start
Tr = Tr[Tr$Species %in% VecOfMammalianSpecies,] 

AGG1 = aggregate(Tr$PercentMatches, by = list(Tr$Species), FUN = mean)
AGG2 = aggregate(Tr$ATCont, by = list(Tr$Species), FUN = mean)
AGG = merge(AGG1,AGG2, by = "Group.1"); names(AGG) = c('Species','MeanPercentMatches','MeanAtCont')
cor.test(AGG$MeanPercentMatches,AGG$MeanAtCont, method = 'spearman') # 
plot(AGG$MeanPercentMatches,AGG$MeanAtCont)

######  ANALYSIS TO THINK FOR FUTURE

# between species analysis:
AGG1 = aggregate(Tr$CopyNumber, by = list(Tr$Species), FUN = mean)
AGG2 = aggregate(Tr$PercentMatches, by = list(Tr$Species), FUN = mean)
AGG = merge(AGG1,AGG2, by = 'Group.1')
names(AGG) = c('Species','CopyNumber','PercentMatches')
plot(AGG$CopyNumber,AGG$PercentMatches)
dev.off()

