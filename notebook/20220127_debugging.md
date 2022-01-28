
# Debugging 20220117 sim run

See last post for list of sims that need debugging:

- 1231717 this simulation supposedly failed at the mutation plotting step
   -   I think this ran out of memory, but I'm not sure why. 21,000 variants. "system call failed: Cannot allocate memory"
   
- 1231101 this simulation supposedly failed at the pca step
  - 44485 SNP variants. "system call failed: Cannot allocate memory"
  - "highly-polygenic_1-trait__SS-Clines_N-equal_m_breaks\n1231101"
  -  


- 1231103 this simulation supposedly failed at the output step
  - this simulation only produced 13000 variants for some reason! 
  - ``` Warning message:
In max(which(a > b * 1.5)) :
  no non-missing arguments to max; returning -Inf
Error in seq_len(min(n, nu)) :
  argument must be coercible to non-negative integer
Calls: lfmm2 -> svd -> La.svd
Execution halted
Warning message:
system call failed: Cannot allocate memory
```
  - 

Before I rerun R, let's delete former products. This code should remove all output files from the R run.
```
cd /work/lotterhos/MVP-NonClinalAF/sim_output_20220117
rm *.pdf
rm *_R.error
rm *_R.out
rm *Rout*
rm *_genotypes*
```

I also edited `e-run_nonAF_sims_1R-fastruns-20220117.sh` to increase the memory for each R run, and run half the jobs at a time.
I'm not sure why this works (I thought each job runs on a core and the max/core is 5Gb, but when I use `seff` it shows that some jobs use up to 20Gb of memory.
There must be something happening under the hood I don't understand!
```
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --array=2-1000%30
```
