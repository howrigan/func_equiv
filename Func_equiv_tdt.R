

## Reading in TDT results from functional equivalence


## Goals:
## - check transmission rate
## - split by different parts of the genome
## - check against other vcf files?

## NOTE: installed BEDTools2 on my computer at: /Users/daniel/Dropbox (Partners HealthCare)/tools/bedtools2
...might still need to get this in the ./bin to run properly

## had to install data.table opening R.app
library(data.table)


## ===================== CHAPTERS:

## 0 - Key HAIL commands
## 0 - first 20K lines testing

## 1 - Singleton Transmission in parent_and_kids.vcf

## 2 - Broad pedigree file

## 3 - singleton count testing

## 4 - de novo finder

## 5 - combine singleton trans, mendelian errors, and de novos

## ===================== END of CHAPTERS:



## 0 - Key HAIL commands

## running TDT
(hc.read('parent_and_kids.vds').split_multi().tdt("vcf/parent_and_kids.fam").export_variants("output/tdt_results.tsv", "Variant = v, va.tdt.*"))


## pulling out singleton lists
vds.filter_variants_expr('gs.map(g => g.nNonRefAlleles()).sum() == 1')
.export_genotypes('[file].txt.bgz', 'variant = v, individual = s')






## 0 - first 20K lines testing


setwd('/Users/daniel/Documents/functional_equivalence/vcf')

zless parent_and_kids.vcf.gz | head -n20000 >  parent_and_kids_first20000lines.vcf
zless parent1_kids2.vcf.gz | head -n20000 >  parent1_kids2_first20000lines.vcf
zless parent2_kids1.vcf.gz | head -n20000 >  parent2_kids1_first20000lines.vcf






## 1 - Singleton Transmission in parent_and_kids.vcf

setwd('/Users/daniel/Documents/functional_equivalence')

tdt <- fread('output/tdt_results.tsv',stringsAsFactors=F)
dim(tdt) ## 24096722


## looks like the chi2 value is incorrect

## Singleton TDT counts:

tdt2 <- subset(tdt,(tdt$nTransmitted==1 & tdt$nUntransmitted==0) | (tdt$nTransmitted==0 & tdt$nUntransmitted==1))
dim(tdt2) ## 3291224

sum(tdt2$nTransmitted) ## 1519306
sum(tdt2$nUntransmitted) ## 1771918

sum(tdt2$nTransmitted) / sum(tdt2$nUntransmitted) ## 0.8574358
sum(tdt2$nTransmitted) / (sum(tdt2$nTransmitted) + sum(tdt2$nUntransmitted)) ## 0.4616234

## create .bed file

tmp <- strsplit(tdt2$Variant,':')
tmp2 <- sapply(tmp,'[',1)
CHROM <- substr(tmp2,start=4,stop=nchar(tmp2))
POS <- sapply(tmp,'[',2)
REF <- sapply(tmp,'[',3)
ALT <- sapply(tmp,'[',4)

tdt2$CHROM <- CHROM
tdt2$POS <- POS
tdt2$REF <- REF
tdt2$ALT <- ALT

## write to file
write.table(tdt2,'output/parent_and_kids_parentalSingletonTransmission.tsv',col=T,row=F,quo=F,sep='\t')





## 2 - Broad pedigree file

setwd('/Users/daniel/Documents/functional_equivalence/vcf')

ped <- read.table('Broad.ped',stringsAsFactors=F)

## recode FID
numb <- as.numeric(as.factor(ped$V3))
ped$V1 <- paste0('fam',numb)

names(ped) <- c('FID','IID','PID','MID','SEX','AFF')

## duplicated IDs?
sum(duplicated(ped$IID)) ## 0
sum(duplicated(ped$FID)) ## 0


father <- cbind.data.frame(ped$FID,ped$PID)
father$PID <- '0'
father$MID <- '0'
father$SEX <- '1'
father$AFF <- '-9'
names(father) <- c('FID','IID','PID','MID','SEX','AFF')
## remove duplicate samples
father2 <- father[!duplicated(father$IID),]


mother <- cbind.data.frame(ped$FID,ped$MID)
mother$PID <- '0'
mother$MID <- '0'
mother$SEX <- '2'
mother$AFF <- '-9'
names(mother) <- c('FID','IID','PID','MID','SEX','AFF')
## remove duplicate samples
mother2 <- mother[!duplicated(mother$IID),]


ped2 <- rbind(ped,mother2,father2)
ped2 <- ped2[order(ped2$FID),]

sum(!duplicated(ped2$IID)) ## 100

table(table(ped2$FID))
 # 3  4 
 # 8 19 

## 8 trios
## 19 quads




## read in sample list from VCFs


## ------ parent_and_kids

## zless parent_and_kids_first5000lines.vcf.bgz | grep '#CHROM' > parent_and_kids.id
## perl -pe 's{\t}{\n}g' parent_and_kids.id > parent_and_kids.id.txt
## remove top lines

pk <- read.table('parent_and_kids.id.txt',stringsAsFactors=F)
dim(pk) ## 96

ped_pk <- subset(ped2,ped2$IID %in% pk$V1)

miss <- subset(ped2,!ped2$IID %in% pk$V1)
miss
#                                               FID                                            IID                                            PID
# 7  Broad__Coriell_Trios_2__H_IJ-HG02818-HG02818_1 Broad__Coriell_Trios_2__H_IJ-HG02818-HG02818_1 Broad__Coriell_Trios_2__H_IJ-HG02816-HG02816_1
# 71 Broad__Coriell_Trios_2__H_IJ-HG02818-HG02818_1 Broad__Coriell_Trios_2__H_IJ-HG02817-HG02817_1                                              0
# 72 Broad__Coriell_Trios_2__H_IJ-HG02818-HG02818_1 Broad__Coriell_Trios_2__H_IJ-HG02816-HG02816_1                                              0
# 13                    Broad__SSC-Batch1__SSC00518                    Broad__SSC-Batch1__SSC00518                    Broad__SSC-Batch1__SSC00523
#                                               MID SEX AFF
# 7  Broad__Coriell_Trios_2__H_IJ-HG02817-HG02817_1   2  -9
# 71                                              0   2  -9
# 72                                              0   1  -9
# 13                    Broad__SSC-Batch1__SSC00522   1  -9

## write to file
write.table(ped_pk,'parent_and_kids.fam',col=F,row=F,quo=F,sep='\t')


## ------ parent1_kids2

## less parent1_kids2_first20000lines.vcf | grep '#CHROM' > parent1_kids2.id
## perl -pe 's{\t}{\n}g' parent1_kids2.id > parent1_kids2.id.txt
## remove top lines

pk1 <- read.table('parent1_kids2.id.txt',stringsAsFactors=F)
dim(pk1) ## 44

ped_pk1 <- subset(ped2,ped2$IID %in% pk1$V1)

## how many families are NOT split?
fam <- names(table(ped_pk1$FID))
parent <- NA
kid <- NA

for (i in 1:length(fam)) {
	pp <- subset(ped_pk1,ped_pk1$FID==fam[i])
	parent[i] <- sum(pp$PID==0)
	kid[i] <- sum(pp$PID!=0)
}

(chk <- cbind.data.frame(fam,parent,kid))
(chk2 <- chk[chk$parent > 0 & chk$kid > 0,])
## no full families

## square up discrepnacies with main ped file
# chk3 <- ped2[ped2$FID %in% chk2$fam,]
# subset(chk3,chk3$IID %in% pk1$V1)

## write to file
write.table(ped_pk1,'parent1_kids2.fam',col=F,row=F,quo=F,sep='\t')


## ------ parent2_kids1

# zless parent2_kids1_first20000lines.vcf | grep '#CHROM' > parent2_kids1.id
# perl -pe 's{\t}{\n}g' parent2_kids1.id > parent2_kids1.id.txt
# remove top lines

pk2 <- read.table('parent2_kids1.id.txt',stringsAsFactors=F)
dim(pk2) ## 52

ped_pk2 <- subset(ped2,ped2$IID %in% pk2$V1)

## how many families are NOT split?
fam <- names(table(ped_pk2$FID))
parent <- NA
kid <- NA

for (i in 1:length(fam)) {
	pp <- subset(ped_pk2,ped_pk2$FID==fam[i])
	parent[i] <- sum(pp$PID==0)
	kid[i] <- sum(pp$PID!=0)
}

(chk <- cbind.data.frame(fam,parent,kid))
(chk2 <- chk[chk$parent > 0 & chk$kid > 0,])
## no full families

## write to file
write.table(ped_pk2,'parent2_kids1.fam',col=F,row=F,quo=F,sep='\t')








## 3 - singleton count testing


## HOW to restrict to COMBINED singletons
## - issue: could be a singleton in one VCF, but not the other
## - will lead to false under-transmission effect
## POSSIBLE solution:
## - read in parent and singleton tdt results
## - restrict PARENTAL singleton lists to these variants
## - if done to both lists, will restrict to singletons




setwd('/Users/daniel/Documents/functional_equivalence')

## ---- parent1_kids2
library(data.table)

## read in singleton transmission list
s_tran <- fread('output/parent_and_kids_parentalSingletonTransmission.tsv',stringsAsFactors=F)
dim(s_tran)
## 3291224

## ---- parent1_kids2
pk1 <- fread('gzip -dc output/parent1_kids2_singletons.txt.bgz',stringsAsFactors=F)
nrow(pk1) ## 6622489

sum(pk1$variant %in% s_tran$Variant) ## 1474025 overlap

pk1_fam <- read.table('vcf/parent1_kids2.fam',stringsAsFactors=F)
nrow(pk1_fam) ## 44

## double check ID matches
ids <- unique(pk1$individual)
sum(ids %in% pk1_fam$V2) / length(ids)
sum(pk1_fam$V2 %in% ids) / length(pk1_fam$V2)

## pull out position
tmp <- strsplit(pk1$variant,':')
pk1$CHROM <- sapply(tmp,'[',1)
# pk1$POS <- as.numeric(sapply(tmp,'[',2)) 
# pk1$REF <- sapply(tmp,'[',3) 
# pk1$ALT <- sapply(tmp,'[',4) 

## adjust pk1 to remove autosomes
pk1 <- subset(pk1,pk1$CHROM != 'chrX')
nrow(pk1) ## 6418641


## split parental and kid singletons
pk1_parents <- subset(pk1_fam,pk1_fam$V3=='0')
s_p1 <- as.data.frame(subset(pk1,pk1$individual %in% pk1_parents$V2))
nrow(s_p1) ## 3937894

pk1_kids <- subset(pk1_fam,pk1_fam$V3!='0')
s_k2 <- as.data.frame(subset(pk1,pk1$individual %in% pk1_kids$V2))
nrow(s_k2) ## 2480747


## subset out variants that aren't parental singletons 
sum(s_p1$variant %in% s_tran$Variant) ## 21207
sum(s_k2$variant %in% s_tran$Variant) ## 1367193
s_96_p1 <- subset(s_p1,s_p1$variant %in% s_tran$Variant)
s_96_k2 <- subset(s_k2,s_k2$variant %in% s_tran$Variant)


chk <- subset(s_p1,s_p1$variant %in% s_tran$Variant)
table(chk$CHROM)
chk <- subset(s_k2,s_k2$variant %in% s_tran$Variant)
table(chk$CHROM)

## Most parental matches coming from the X chromosome




## ---- parent2_kids1
pk2 <- fread('gzip -dc output/parent2_kids1_singletons.txt.bgz',stringsAsFactors=F)
pk2_fam <- read.table('vcf/parent2_kids1.fam',stringsAsFactors=F)

## double check ID matches
ids <- unique(pk2$individual)
sum(ids %in% pk2_fam$V2) / length(ids)
sum(pk2_fam$V2 %in% ids) / length(pk2_fam$V2)

## pull out position
tmp <- strsplit(pk2$variant,':')
pk2$CHROM <- sapply(tmp,'[',1)
# pk2$POS <- as.numeric(sapply(tmp,'[',2)) 
# pk2$REF <- sapply(tmp,'[',3) 
# pk2$ALT <- sapply(tmp,'[',4) 


## adjust pk2 to remove autosomes
pk2 <- subset(pk2,pk2$CHROM != 'chrX')
nrow(pk2) ## 6354637


## split parental and kid singletons
pk2_parents <- subset(pk2_fam,pk2_fam$V3=='0')
s_p2 <- as.data.frame(subset(pk2,pk2$individual %in% pk2_parents$V2))
nrow(s_p2) ## 4462535

pk2_kids <- subset(pk2_fam,pk2_fam$V3!='0')
s_k1 <- as.data.frame(subset(pk2,pk2$individual %in% pk2_kids$V2))
nrow(s_k1) ## 1892102


## subset out variants that aren't parental singletons 
sum(s_p2$variant %in% s_tran$Variant) ## 2827205
sum(s_k1$variant %in% s_tran$Variant) ## 4641
s_96_p2 <- subset(s_p2,s_p2$variant %in% s_tran$Variant)
s_96_k1 <- subset(s_k1,s_k1$variant %in% s_tran$Variant)



chk <- subset(s_p2,s_p2$variant %in% s_tran$Variant)
table(chk$CHROM)
chk <- subset(s_k1,s_k1$variant %in% s_tran$Variant)
table(chk$CHROM)



## VERDICT: 
 - Parent and Kids 1 have very few singletons
 - Parent and Kids 2 have a lot more




## Loop through kids2
 # - grab singleton count
 # - find child
 # - match variant 
 # - grab transmission count  



## ------ Kids1 Blocking scheme

pk2_kids <- subset(pk2_fam,pk2_fam$V3!='0')
names(pk2_kids) <- c('FID','IID','PID','MID','SEX','AFF')

father_shet_count <- NA
mother_shet_count <- NA
father_matched <- NA
mother_matched <- NA
father_tr <- NA
mother_tr <- NA
parental_tr <- NA

for (i in 1:nrow(pk2_kids)) {

	# klist <- subset(s_96_k1,s_96_k1$individual == pk2_kids$IID[i])
	# flist <- subset(s_96_p1,s_96_p1$individual == pk2_kids$PID[i])
	# mlist <- subset(s_96_p1,s_96_p1$individual == pk2_kids$MID[i])

	klist <- subset(s_k1,s_k1$individual == pk2_kids$IID[i])
	flist <- subset(s_p1,s_p1$individual == pk2_kids$PID[i])
	mlist <- subset(s_p1,s_p1$individual == pk2_kids$MID[i])

	## get parental counts
	father_shet_count[i] <- nrow(flist)
	mother_shet_count[i] <- nrow(mlist)
	## get tranmission count
	father_matched[i] <- sum(klist$variant %in% flist$variant)
	mother_matched[i] <- sum(klist$variant %in% mlist$variant)
	## get transmission rate
	father_tr[i] <- father_matched[i] / father_shet_count[i]
	mother_tr[i] <- mother_matched[i] / mother_shet_count[i]
	parental_tr[i] <- (father_matched[i] + mother_matched[i]) / (father_shet_count[i] + mother_shet_count[i])

	print(i)

} ## END of i LooP

dat <- cbind.data.frame(pk2_kids,
father_shet_count,
mother_shet_count,
father_matched,
mother_matched,
father_tr,
mother_tr,
parental_tr)

dat[,c('IID','PID','MID','father_matched','mother_matched','parental_tr')]

## write to file
write.table(dat,'output/parent1_kids1_singletonTransmission_summary.txt',col=T,row=F,quo=F,sep='\t')





## ------ Kids2 Blocking scheme

pk1_kids <- subset(pk1_fam,pk1_fam$V3!='0')
names(pk1_kids) <- c('FID','IID','PID','MID','SEX','AFF')

father_shet_count <- NA
mother_shet_count <- NA
father_matched <- NA
mother_matched <- NA
father_tr <- NA
mother_tr <- NA
parental_tr <- NA

for (i in 1:nrow(pk2_kids)) {

	# klist <- subset(s_96_k2,s_96_k2$individual == pk1_kids$IID[i])
	# flist <- subset(s_96_p2,s_96_p2$individual == pk1_kids$PID[i])
	# mlist <- subset(s_96_p2,s_96_p2$individual == pk1_kids$MID[i])

	klist <- subset(s_k2,s_k2$individual == pk1_kids$IID[i])
	flist <- subset(s_p2,s_p2$individual == pk1_kids$PID[i])
	mlist <- subset(s_p2,s_p2$individual == pk1_kids$MID[i])

	## get parental counts
	father_shet_count[i] <- nrow(flist)
	mother_shet_count[i] <- nrow(mlist)
	## get tranmission count
	father_matched[i] <- sum(klist$variant %in% flist$variant)
	mother_matched[i] <- sum(klist$variant %in% mlist$variant)
	## get transmission rate
	father_tr[i] <- father_matched[i] / father_shet_count[i]
	mother_tr[i] <- mother_matched[i] / mother_shet_count[i]
	parental_tr[i] <- (father_matched[i] + mother_matched[i]) / (father_shet_count[i] + mother_shet_count[i])

print(i)

} ## END of i LooP

dat <- cbind.data.frame(pk1_kids,
father_shet_count,
mother_shet_count,
father_matched,
mother_matched,
father_tr,
mother_tr,
parental_tr)

dat[,c('IID','PID','MID','father_matched','mother_matched','parental_tr')]


## write to file
write.table(dat,'output/parent2_kids2_singletonTransmission_summary.txt',col=T,row=F,quo=F,sep='\t')







## 4 - de novo finder

## Run de_novo_finder.py

python /Users/daniel/Documents/functional_equivalence/de_novo_finder_3.py \
/Users/daniel/Documents/functional_equivalence/vcf/parent_and_kids.vcf.gz \
/Users/daniel/Documents/functional_equivalence/vcf/parent_and_kids.fam \
/Users/daniel/Dropbox/denovo/all_ESP_counts_5.28.13.txt \
> /Users/daniel/Documents/functional_equivalence/output/parent_and_kids.deNovo





## 5 - combine singleton trans, mendelian errors, and de novos

setwd('/Users/daniel/Documents/functional_equivalence/output')

p1 <- read.table('parent1_kids1_singletonTransmission_summary.txt',h=T,stringsAsFactors=F)
p1$scheme <- '1'
p2 <- read.table('parent2_kids2_singletonTransmission_summary.txt',h=T,stringsAsFactors=F)
p2$scheme <- '2'

trans <- rbind.data.frame(p1,p2)

## Mendelian errors
imend <- read.table('parent_and_kids.imendel',h=T,stringsAsFactors=F)
names(imend) <- c('FID','IID','N_MENDEL_ERRORS','NSNP_MENDEL_ERRORS')

tmp <- merge(trans,imend[,c('IID','N_MENDEL_ERRORS','NSNP_MENDEL_ERRORS')],by='IID',all.x=T)

## denovos
dnm <- read.table('parent_and_kids.deNovo',h=T,stringsAsFactors=F)

tbl <- table(dnm$Child_ID)
dnm_count <- cbind.data.frame(names(tbl),as.numeric(tbl))
names(dnm_count) <- c('IID','denovo_count_unfiltered')

tmp2 <- merge(tmp,dnm_count,by='IID',all.x=T)

## write to file
write.table(tmp2,'family_transmission_denovo_summary.txt',col=T,row=F,quo=F,sep='\t')



## ---- Quick analysis

cor(tmp2[,c('parental_tr','N_MENDEL_ERRORS','NSNP_MENDEL_ERRORS','denovo_count_unfiltered')])
#                         parental_tr N_MENDEL_ERRORS NSNP_MENDEL_ERRORS denovo_count_unfiltered
# parental_tr               1.0000000      -0.4074415         -0.5200937               0.7732809
# N_MENDEL_ERRORS          -0.4074415       1.0000000          0.9764885              -0.2757560
# NSNP_MENDEL_ERRORS       -0.5200937       0.9764885          1.0000000              -0.3570789
# denovo_count_unfiltered   0.7732809      -0.2757560         -0.3570789               1.0000000

t.test(tmp2$parental_tr[tmp2$scheme==1],tmp2$parental_tr[tmp2$scheme==2])
# t = -3.8716, df = 23.363, p-value = 0.0007565
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.13223076 -0.04018452
# sample estimates:
# mean of x mean of y 
# 0.2089896 0.2951973 



## ---- check within SSC
aa <- strsplit(tmp2$IID,'-')
bb <- sapply(aa,'[',1)

ssc <- subset(tmp2,bb=='Broad__SSC')

t.test(ssc$parental_tr[ssc$scheme==1],ssc$parental_tr[ssc$scheme==2])
# t = -1.6853, df = 17.555, p-value = 0.1096
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.058376817  0.006461236
# sample estimates:
# mean of x mean of y 
# 0.2089896 0.2349474 










































