rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

user = 'Alina'
if(user == 'Alina'){
  setwd('/home/aragret/Alina/COMPARATIVE_GENOMICS')
}

tr = read.table('Derived_data/TRFinder.txt', sep='\t', header=TRUE)

hist(tr$PercentMatches)

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",as.character(s))
  return (nchar(s) - nchar(s2))
}

for(i in 1:nrow(tr)){
  # i = 1
  tr$A[i] = countCharOccurrences('A', as.character(tr$RepeatsRegion[i]))
  tr$G[i] = countCharOccurrences('G', as.character(tr$RepeatsRegion[i]))
  tr$C[i] = countCharOccurrences('C', as.character(tr$RepeatsRegion[i]))
  tr$Ti[i] = countCharOccurrences('T', as.character(tr$RepeatsRegion[i]))
  tr$ATCont[i] = (tr$A[i] + tr$Ti[i]) / nchar(as.character(tr$RepeatsRegion[i]))
  tr$GCCont[i] = (tr$G[i] + tr$C[i]) / nchar(as.character(tr$RepeatsRegion[i]))
}

young_repeats = tr[tr$PercentMatches >= 95,]
old_repeats = tr[tr$PercentMatches < 95,]

# write.table(young_repeats, 'Derived_data/YoungTandemRepeats.txt', sep='\t', row.names = F)
# write.table(old_repeats, 'Derived_data/OldTandemRepeats.txt', sep='\t', row.names = F)


library(ggplot2)
library(gridExtra)

ggplot(tr, aes(ATCont)) + geom_histogram(col="red", fill="green", alpha = .2)

plot1 = ggplot(young_repeats, aes(ATCont)) + geom_histogram(col="red", fill="green", alpha = .2)
plot2 = ggplot(old_repeats, aes(ATCont)) + geom_histogram(col="red", fill="green", alpha = .2)
grid.arrange(plot1, plot2, ncol=2)

plot1 = ggplot(young_repeats, aes(GCCont)) + geom_histogram(col="red", fill="green", alpha = .2)
plot2 = ggplot(old_repeats, aes(GCCont)) + geom_histogram(col="red", fill="green", alpha = .2)
grid.arrange(plot1, plot2, ncol=2)

wilcox.test(young_repeats$ATCont, old_repeats$ATCont)
t.test(young_repeats$ATCont, old_repeats$ATCont)

plot1 = ggplot(young_repeats, aes('1', ATCont)) + geom_boxplot(col="red", fill="green", alpha = .2)
plot2 = ggplot(old_repeats, aes('2', ATCont)) + geom_boxplot(col="red", fill="green", alpha = .2)
grid.arrange(plot1, plot2, ncol=2)

##############################################################################
### repeats within each species

length(unique(tr$Species)) # 2040
length(unique(young_repeats$Species)) # 1311
length(unique(old_repeats$Species)) # 1295

one_line = c()
for(sp in unique(tr$Species)){
  young_sp = young_repeats[young_repeats$Species == as.character(sp),]
  old_sp = old_repeats[old_repeats$Species == as.character(sp),]
  if((nrow(old_sp) > 1) & (nrow(young_sp) > 1)){
    a = merge(young_sp, old_sp, by='Species')
    a = a[, c('ATCont.x', 'ATCont.y')]
    test = wilcox.test(a[,1], a[,2])
    one_line = rbind(one_line, c(sp, nrow(young_sp), nrow(old_sp), test$statistic, test$p.value))
  }
}

result = as.data.frame(one_line)
names(result) = c('Species', 'YoungRepeatsNumber', 'OldSpeciesNumber', 'W', 'PValue')
write.table(result, 'Derived_data/ATContentInYoungOldTR.txt', sep='\t', quote=FALSE,
            row.names = FALSE)

###########################################################################
### 