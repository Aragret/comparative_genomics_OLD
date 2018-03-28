setwd('/home/aragret/Alina/COMPARATIVE_GENOMICS')

repeats = read.table('Derived_data/AllRepeatsBaseContent.txt', sep='\t', header=TRUE)
CHOR = read.table('Derived_data/MitGenomics.txt', sep='\t', header=TRUE)

CHOR = CHOR[, c('Species', 'TAXON', 'GCCont', 'ECO.Maximum.longevity..yrs.',
                'ECO.Female.maturity..days.', 'GenomeLength', 
                'REP.LengthOfTandemRepeats', 'REP.LengthOfTandemRepeatsWithoutOverlaps')]

for (i in 1:nrow(repeats)){
  repeats$dir_GC[i] = (repeats$dir_G[i] + repeats$dir_C[i]) / 
    (repeats$dir_A[i] + repeats$dir_T[i] + repeats$dir_G[i] + repeats$dir_C[i])
  repeats$symm_GC[i] = (repeats$symm_G[i] + repeats$symm_C[i]) / 
    (repeats$symm_A[i] + repeats$symm_T[i] + repeats$symm_G[i] + repeats$symm_C[i])
  repeats$compl_GC[i] = (repeats$compl_G[i] + repeats$compl_C[i]) / 
    (repeats$compl_A[i] + repeats$compl_T[i] + repeats$compl_G[i] + repeats$compl_C[i])
  repeats$inv_GC[i] = (repeats$inv_G[i] + repeats$inv_C[i]) / 
    (repeats$inv_A[i] + repeats$inv_T[i] + repeats$inv_G[i] + repeats$inv_C[i])
  repeats$tand_GC[i] = (repeats$tand_G[i] + repeats$tand_C[i]) / 
    (repeats$tand_A[i] + repeats$tand_T[i] + repeats$tand_G[i] + repeats$tand_C[i])
  
}

df = merge(CHOR, repeats[, c('Species', 'dir_GC', 'symm_GC', 'compl_GC', 'inv_GC', 'tand_GC')], 
           by='Species')

for (k in unique(df$TAXON)){
  
  TempData = df[df$TAXON == k,]
  
  TempData$dirLog2Ration = log2(TempData$GCCont/TempData$dir_GC)
  Median = median(TempData$dirLog2Ration)
  pvalue = as.numeric(wilcox.test(TempData$GCCont, TempData$dir_GC, paired = TRUE)[3])
  print(c(k, Median, pvalue))
  # median - [1] 3.11528
  # p-value 0
  
  TempData$symmLog2Ration = log2(TempData$GCCont/TempData$symm_GC)
  Median = median(TempData$symmLog2Ration)
  pvalue = as.numeric(wilcox.test(TempData$GCCont, TempData$symm_GC, paired = TRUE)[3])
  # median - [1] 3.187956
  # p-value 0
  print(c(k, Median, pvalue))
  
  TempData$complLog2Ration = log2(TempData$GCCont/TempData$compl_GC)
  Median = median(TempData$complLog2Ration)
  pvalue = as.numeric(wilcox.test(TempData$GCCont, TempData$compl_GC, paired = TRUE)[3])
  # median - [1] 3.233705
  # p-value 0
  print(c(k, Median, pvalue))
  
  TempData$invLog2Ration = log2(TempData$GCCont/TempData$inv_GC)
  Median = median(TempData$invLog2Ration)
  pvalue = as.numeric(wilcox.test(TempData$GCCont, TempData$inv_GC, paired = TRUE)[3])
  # median - [1] 3.103899
  # p-value 0
  print(c(k, Median, pvalue))
  
  TempData$tandLog2Ration = log2(TempData$GCCont/TempData$tand_GC)
  Median = median(TempData$tandLog2Ration[!is.na(TempData$tandLog2Ration)]) 
  pvalue = as.numeric(wilcox.test(TempData$GCCont, TempData$tand_GC, paired = TRUE)[3])
  # median - [1] 0.5582295
  # p-value < 2.233208e-157
  print(c(k, Median, pvalue))
  
}
df$dirLog2Ration = log2(df$GCCont/df$dir_GC)
Median = median(df$dirLog2Ration);  Median 
pvalue = as.numeric(wilcox.test(df$GCCont, df$dir_GC, paired = TRUE)[3]); pvalue
# median - [1] 3.11528
# p-value 0

df$symmLog2Ration = log2(df$GCCont/df$symm_GC)
Median = median(df$symmLog2Ration);  Median 
pvalue = as.numeric(wilcox.test(df$GCCont, df$symm_GC, paired = TRUE)[3]); pvalue
# median - [1] 3.187956
# p-value 0

df$complLog2Ration = log2(df$GCCont/df$compl_GC)
Median = median(df$complLog2Ration);  Median 
pvalue = as.numeric(wilcox.test(df$GCCont, df$compl_GC, paired = TRUE)[3]); pvalue
# median - [1] 3.233705
# p-value 0

df$invLog2Ration = log2(df$GCCont/df$inv_GC)
Median = median(df$invLog2Ration);  Median 
pvalue = as.numeric(wilcox.test(df$GCCont, df$inv_GC, paired = TRUE)[3]); pvalue
# median - [1] 3.103899
# p-value 0

df$tandLog2Ration = log2(df$GCCont/df$tand_GC)
Median = median(df$tandLog2Ration[!is.na(df$tandLog2Ration)]);  Median 
pvalue = as.numeric(wilcox.test(df$GCCont, df$tand_GC, paired = TRUE)[3]); pvalue
# median - [1] 0.5582295
# p-value < 2.233208e-157

####
TempData = df[df$TAXON == 'Mammalia',]
summary(TempData)

med = median(TempData$ECO.Maximum.longevity..yrs.[!is.na(TempData$ECO.Maximum.longevity..yrs.)])
long_lived = TempData[TempData$ECO.Maximum.longevity..yrs. > med,]
short_lived = TempData[TempData$ECO.Maximum.longevity..yrs. <= med,]

long_lived$tandLog2Ration = log2(long_lived$GCCont/long_lived$tand_GC)
Median = median(long_lived$tandLog2Ration[!is.na(long_lived$tandLog2Ration)]);  Median 
# [1] -0.03713389
pvalue = as.numeric(wilcox.test(long_lived$GCCont, long_lived$tand_GC, paired = TRUE)[3]); pvalue
# [1] 0.7984496

short_lived$tandLog2Ration = log2(short_lived$GCCont/short_lived$tand_GC)
Median = median(short_lived$tandLog2Ration[!is.na(short_lived$tandLog2Ration)]);  Median 
# -0.02661405
pvalue = as.numeric(wilcox.test(short_lived$GCCont, short_lived$tand_GC, paired = TRUE)[3]); pvalue
# [1] 0.7534716

quantile_vec = quantile(TempData$ECO.Maximum.longevity..yrs., probs = c(0.25, 0.5, 0.75),
                        na.rm = TRUE)

long_lived = TempData[TempData$ECO.Maximum.longevity..yrs. > quantile_vec[3],]
short_lived = TempData[TempData$ECO.Maximum.longevity..yrs. <= quantile_vec[1],]

long_lived$tandLog2Ration = log2(long_lived$GCCont/long_lived$tand_GC)
Median = median(long_lived$tandLog2Ration[!is.na(long_lived$tandLog2Ration)]);  Median 
# -0.2226203
pvalue = as.numeric(wilcox.test(long_lived$GCCont, long_lived$tand_GC, paired = TRUE)[3]); pvalue
# 0.03874834

short_lived$tandLog2Ration = log2(short_lived$GCCont/short_lived$tand_GC)
Median = median(short_lived$tandLog2Ration[!is.na(short_lived$tandLog2Ration)]);  Median 
# -0.03793129
pvalue = as.numeric(wilcox.test(short_lived$GCCont, short_lived$tand_GC, paired = TRUE)[3]); pvalue
# 0.8309879

almost_short_lived = TempData[TempData$ECO.Maximum.longevity..yrs. > quantile_vec[1] &
                                TempData$ECO.Maximum.longevity..yrs. < quantile_vec[2],]

almost_long_lived = TempData[TempData$ECO.Maximum.longevity..yrs. >= quantile_vec[2] &
                                TempData$ECO.Maximum.longevity..yrs. < quantile_vec[3],]

almost_long_lived$tandLog2Ration = log2(almost_long_lived$GCCont/almost_long_lived$tand_GC)
Median = median(almost_long_lived$tandLog2Ration[!is.na(almost_long_lived$tandLog2Ration)]);  Median 
# 0.03842235
pvalue = as.numeric(wilcox.test(almost_long_lived$GCCont, almost_long_lived$tand_GC, paired = TRUE)[3]); pvalue
# 0.07350794

almost_short_lived$tandLog2Ration = log2(almost_short_lived$GCCont/almost_short_lived$tand_GC)
Median = median(almost_short_lived$tandLog2Ration[!is.na(almost_short_lived$tandLog2Ration)]);  Median 
# -0.03090958
pvalue = as.numeric(wilcox.test(almost_short_lived$GCCont, almost_short_lived$tand_GC, paired = TRUE)[3]); pvalue
# 0.9474012

#########################################################################
## TR in mammals

CHOR$TRCoverage = CHOR$REP.LengthOfTandemRepeats / CHOR$GenomeLength
summary(CHOR$TRCoverage)

summary(CHOR$REP.LengthOfTandemRepeatsWithoutOverlaps / CHOR$GenomeLength)

mamm = CHOR[CHOR$TAXON == 'Mammalia',]

temp_data = mamm[!is.na(mamm$ECO.Maximum.longevity..yrs.),]
cor.test(temp_data$TRCoverage, temp_data$ECO.Maximum.longevity..yrs., method = 'spearman')

for(k in taxon_vec){
  TempData = CHOR[CHOR$TAXON == k,]
  print(k)
  print(cor.test(log(TempData$TRCoverage, base = exp(1)), log(TempData$ECO.Female.maturity..days., base = exp(1)), method = 'spearman')
)
}

by(CHOR, CHOR$TAXON, nrow)
table(CHOR$TAXON)

### nice! Check it with PIC

features = c('TRCoverage', 'ECO.Female.maturity..days.')
row.names(CHOR) = CHOR$Species
CHOR$TRCoverage = CHOR$REP.LengthOfTandemRepeats / CHOR$GenomeLength
summary(CHOR$TRCoverage)

one_line = c()
for(k in taxon_vec){
  # k = 'Actinopterygii'
  TempData = data[data$TAXON == k, features]
  TempData = TempData[!is.na(TempData[,1]),]; TempData = TempData[!is.na(TempData[,2]),]
  temp_diff = setdiff(tree2$tip.label, row.names(TempData))
  temp_tree <- drop.tip(tree2, temp_diff)
  contrasts <- apply(TempData, 2, pic, temp_tree)
  a = cor.test(scale(contrasts[,1]), scale(contrasts[,2]), method = 'spearman')
  # print(k)
  one_line = rbind(one_line, c(k, a$estimate, a$p.value))
}
one_line
