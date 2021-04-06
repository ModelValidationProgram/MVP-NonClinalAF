# 20210315


* The VCF files are too large to deal with, so I subsampled 200 random individuals (which should give us on average ~ 100 per population) and an MAF of 0.01

`vcftools --vcf 3384725_InversionVCF.vcf --max-indv 200 --maf 0.05 --out 3384725_InversionVCF_200ind_MAF05 --recode --recode-INFO-all`

The `--recode --recode-INFO-all` tells it to output a vcf file AND keep all the information.

Sara and I need to discuss the MAF filter - we've talked about 0.01 but are using 0.05 for now

```
Parameters as interpreted:
	--vcf 3384725_InversionVCF.vcf
	--recode-INFO-all
	--maf 0.05
	--max-indv 200
	--out 3384725_InversionVCF_200ind_MAF05
	--recode