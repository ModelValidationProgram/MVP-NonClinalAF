# Vcf file conversion check


Katie wasn't sure this line of code:

```
sed 's/\.		/\.	0	/g' ${outpath}${i}"_plusneut_MAF01.recode.vcf" > ${outpath}${i}"_plusneut_MAF01.recode2.vcf" 
```

was actually replacing a space with a dot in the vcf file produced by pyslim. I confirmed that it was behaving properly by reading it into R using `vcfR`. 

Screenshot:

![](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/img/Screenshot%20from%202021-06-30%2016-11-55.png)
