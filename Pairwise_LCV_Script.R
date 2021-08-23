#!/usr/bin/env Rscript

library(data.table)

setwd('/data/eg03/rjabalamel/Mendelian_Randomization/summary_stats/medication_use')

ldf <- list() #Create an empty list to append all summary statistics in the PWD.
list_sumstat <- list.files(pattern = "*sumstats.gz")

target_comb <- data.frame(t(combn(list_sumstat, 2))) # Creat a DF for all pairwise combination of traits.

LCV_matrix<- matrix(nrow= nrow(target_comb), ncol= 12, data= NA)

for (i in 1:nrow(target_comb)){

  #Start with data munged using the ldsc package

  trait1File= paste("/data/eg03/rjabalamel/Mendelian_Randomization/summary_stats/medication_use/",
                    target_comb$X1[[i]],sep="")

  trait2File= paste("/data/eg03/rjabalamel/Mendelian_Randomization/summary_stats/medication_use/",
                    target_comb$X2[[i]],sep="")

  #Load trait 1 data and calculate Zs
  d1 = na.omit(read.table(gzfile(trait1File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

  #Load trait 2 data and calculate Zs
  d2 = na.omit(read.table(gzfile(trait2File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

  #Load LD scores
  setwd('/data/eg03/rjabalamel/LCV/EU_LDscores/eur_w_ld_chr/')
  LDfiles <- list.files(path =".", pattern = "ldscore.gz")
  temp <- lapply(LDfiles, fread, sep="\t")
  LDscores_unsorted <- rbindlist( temp )
  d3= LDscores_unsorted[with(LDscores_unsorted, order(CHR, BP)), ]


  #Merge
  m = merge(d3,d1,by="SNP")
  data = merge(m,d2,by="SNP")

  #Sort by position
  #data = data[order(data[,"CHR"],data[,"BP"]),]
  data = data[with(data, order(CHR, BP)), ]

  #Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
  #If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
  mismatch = which(data$A1.x!=data$A1.y,arr.ind=TRUE)
  data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
  data[mismatch,]$A1.y = data[mismatch,]$A1.x
  data[mismatch,]$A2.y = data[mismatch,]$A2.x


  #Run LCV-need to setwd to directory containing LCV package
  setwd("/data/eg03/rjabalamel/LCV/R")
  source("/data/eg03/rjabalamel/LCV/R/RunLCV.R")

  try({LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
    LCV_matrix[i,]<- c(target_comb$X1[[i]], target_comb$X2[[i]], LCV$zscore, LCV$pval.gcpzero.2tailed,
              LCV$gcp.pm, LCV$gcp.pse, LCV$rho.est, LCV$rho.err,
              LCV$pval.fullycausal[[1]], LCV$pval.fullycausal[[2]],
              LCV$h2.zscore[[1]], LCV$h2.zscore[[2]])
    })

}

LCV_results <- data.frame(trait_1= LCV_matrix[,1], trait_2= LCV_matrix[,2], zscore= LCV_matrix[,3],
                          gcp_p.Value= LCV_matrix[,4], gcp.pm= LCV_matrix[,5],
                          gcp.pse= LCV_matrix[,6], rho.est= LCV_matrix[,7],
                          rho.err= LCV_matrix[,8], pval.fullycausal_1= LCV_matrix[,9],
                          pval.fullycausal_2= LCV_matrix[,10],
                          h2.zscore_t1= LCV_matrix[,11],
                          h2.zscore_t2= LCV_matrix[,12])

final_DF <- cbind(target_comb, LCV_results)


write.table(final_DF, file = "/data/eg03/rjabalamel/Mendelian_Randomization/LCV_results/medication_use/medication_and_drug_LCV_resuls.txt",
            sep='\t', row.names=FALSE)
                          
rm(list = ls())
