## Ku=ku

rm(list=ls(all=TRUE))
user = 'Kostya'
if (user == 'Kostya')
{
YoungTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/YoungTandemRepeats.txt', header = TRUE, sep = '\t')
OldTr = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/2_DERIVED/OldTandemRepeats.txt', header = TRUE, sep = '\t')
}

YoungAgg = aggregate(YoungTr$ATCont, by = list(YoungTr$Species), FUN = mean)
names(YoungAgg) = c('Species','YoungTrMeanAt')

OldAgg = aggregate(OldTr$ATCont, by = list(OldTr$Species), FUN = mean)
names(OldAgg) = c('Species','OldTrMeanAt')

par(mfrow=c(2,2))
Tr = merge(YoungAgg,OldAgg, by = 'Species')
nrow(Tr) # 566
hist(Tr$YoungTrMeanAt, breaks = 50); 
hist(Tr$OldTrMeanAt, breaks = 50); 

Tr$YoungAtMinusOldAt = Tr$YoungTrMeanAt - Tr$OldTrMeanAt
summary(Tr$YoungAtMinusOldAt)
hist(Tr$YoungAtMinusOldAt, breaks = 100)
wilcox.test(Tr$YoungAtMinusOldAt, mu = 0)
wilcox.test(Tr$YoungTrMeanAt,Tr$OldTrMeanAt, paired = TRUE)
t.test(Tr$YoungTrMeanAt,Tr$OldTrMeanAt, paired = TRUE)

