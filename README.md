# COVID19 MR project
Script for generating the PRS scores & Mendelian Randomisation analysis for the collaborative COVID-19 research.

## NOTE

Please note that the variant IDs in the target BED shoud follow `CHR(numeric):BP(numeric)` [i.e. `19:971949`].
If the target BED contains the rsIDs as variant identifier, update the IDs using PLINK 2.00:

``` plink2 --bfile [target_BED] --snps-only --set-all-var-ids '@:#' --make-bed --out [modified_target_BED] ```


--------------------------------------
## Prerequisite

[PRSice-2](https://www.prsice.info/)

R version 3.2.3 or higher

PLINK 2.00 or higher

## Basic usage

```bash
git clone https://github.com/RezaJF/COVID19_PRS.git

cd COVID19_PRS

./PRS_calculation.sh [PATH_to_PRSice] [PATH_to_imputed_BEDs] [OUTPUT_DIRECTORY]
```
- For accessing accompanying summary statistic data, please contact directly!

## Contact
For more update and instructions email me at: 
reza.jabal@einsteinmed.org or mjaf1d14@soton.ac.uk

