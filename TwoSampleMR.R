library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MRPRESSO)

# Exposure data for all SNPs p < 1e-8
exposure_sumstats <- "/Volumes/BIOBANK/Medication_intake_summary_statistics/psychiatric_genomic_consortium/opioid-exposed_vs._opioid-unexposed_controls_in_European-ancestry_cohorts.txt.gz"
expo_sumstats <- fread(exposure_sumstats)
expo_sumstats$PHEN <- "Opioid_exposure"
exposure_dat <- format_data(expo_sumstats,
                            type="exposure",
                            snp_col = "rsID",
                            beta_col = "Beta",
                            se_col = "SE",
                            effect_allele_col = "Allele1",
                            other_allele_col = "Allele2",
                            #eaf_col = "FRQ",
                            samplesize_col="Total_N",
                            pval_col = "P-value",
                            phenotype_col = "PHEN",
                            min_pval = 1e-8)

# Outcome data
outcome_sumstats <- "/Volumes/BIOBANK/ldsc/Summary_statistics/HGI_COVID-19/COVID19_HGI_A2_ALL_20201020.10k.txt.gz"

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outcome_sumstats,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "all_inv_var_meta_beta",
  se_col = "all_inv_var_meta_sebeta",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "all_meta_AF",
  pval_col = "all_inv_var_meta_p",
  samplesize_col = "all_meta_sample_N",
  min_pval = 1e-8
) %>% 
  mutate(outcome = "severe_covid_A1")


# Harmonize data
dat <- harmonise_data(
  exposure_dat = exposure_dat, 
  outcome_dat = outcome_dat
)

# LD clump
dat_clumped <- clump_data(dat,clump_r2=0.1)

# Perform MR
methods_to_use <- c("mr_egger_regression","mr_ivw","mr_weighted_median")
res <- generate_odds_ratios(mr(dat_clumped, method_list=methods_to_use))
res

# Perform PRESSO
mr_presso(data = dat_clumped,
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
          NbDistribution = 10000,  SignifThreshold = 0.05)

# Sensitivity analysis
# Heterogeneity statistics
mr_heterogeneity(dat_clumped, method_list=c("mr_ivw", "mr_egger_regression"))

# Horizontal pleiotropy
mr_pleiotropy_test(dat_clumped)

# Sinle SNP analysis
res_singleSNP <- mr_singlesnp(dat_clumped, all_method=c("mr_ivw", "mr_egger_regression"))

# Leave-one-out analysis
mr_leaveoneout(dat_clumped, method = mr_ivw)
mr_leaveoneout_plot(mr_leaveoneout(dat_clumped, method = mr_ivw))


# Forest plot
mr_forest_plot(res_singleSNP)

# Funnel plot
mr_funnel_plot(res_singleSNP)

mr_report(
  dat_clumped,
  output_path = "/Volumes/BIOBANK/TwoSampleMR/MRresult",
  output_type = "html",
  author = "Reza Jabal",
  study = "Covid-19 Two Sample MR"
)

rm(list = ls())
