Below is a description of metadata for each set of output tables. The first number is the seed for the simulation.

### 2549986039929_fitnessmat_ind.txt
* an n_deme x n_individual table that indicates the fitness of each individual (in columns) in every metapopulation deme (in rows)

### 2549986039929_fitnessmat.txt
* an n_deme x n_deme table that indicates the mean fitness of individuals from the source deme (in columns) in the transplant deme (in rows) CHECK THIS IS CORRECT

### 2549986039929_ind.txt
Individual information output last generation

* `seed` simulation seed
* `indID` individual ID in SliM. Corresponds to rows in VCF file.
* `indSubpopIndex` index of that individual within the subpopulation
* `subpop` subpopulation ID
* `phen_sal` individual salinity phenotype
* `phen_temp` individual temperature phenotype
* `sal_opt` salinity optimum for that subpopulation
* `temp_opt` temperature optimum for that subpopulation
* `fitness` fitness of the individual in it's subpopulation
    
### 2549986039929_info.txt
Simulation information

* `seed`  simulation seed
* `sim_type` the name of the simulation file used to produce the results
* `MIG_x` migration across longitude populations
* `MIG_y` migration across latitude populations
* `N` population size
* `mu` mutation rate
* `r` recombination rate
* `BURNIN_1` length of first burnin (typically 0 opimums across metapopulation)
* `BURNIN_2` length of second burnin (typically a transition from 0 opimums across metapopulation to final environment values)
* `TOTAL_GEN` the total number of generations
* `METAPOP_SIDE_x` the number of metapopulations across longitude
* `METAPOP_SIDE_y` the number of metapopulations across latitude
* `S` Not used in this simulation.
* `SIGMA_K` Strength of stabilizing selection
* `SIGMA_QTN` Variation in QTN mutation rate

### 2549986039929_LA.txt
This file tracks the amount of local adaptation in the simulation through time

* `seed` simulation seed
* `gen` generation
* `sympatric`  average fitness of demes in sympatry
* `allopatric` average fitness of demes in allopatry
* `local_adapt` amount of local adaptation (`allopatric` - `sympatric` )
* `mean_phen0` mean of first phenotype **WHAT IS THIS**
* `mean_phen1` mean of second phenotype **WHAT IS THIS**
* `cor_sal_popmean` correlation between salinity and population mean phenotype (often inflated due to central limit theorem)
* `cor_temp_popeman` correlation between temperature and population mean phenotype (often inflated due to central limit theorem)
* `cor_sal_ind` correlation between salinity and individual phenotypes
* `cor_temp_ind` correlation between temperature and individual phenotypes

### 2549986039929_mig_mat.txt
This is the migration matrix that was used in the sims. It is used in pyslim for recapitation. Without it, pyslim will run forever and never recapitate.

### 2549986039929_muts.txt
This is a table of information about each mutation that is simulated SLIM

* `seed`  simulation seed
* `mutID` mutation ID - should match mutID in other tables and VCF file
* `muttype` NEED TO ADD TO OUTPUT **WHAT IS THIS**
* `MAF` allele frequency of derived mutation - IS THIS ACTUALLY MAF? **WHAT IS THIS**
* `cor_sal` correlation of mutation allele frequency and salinity at end of simulation (all individuals) - IS THIS CALCULATED BY EACH DEME OR ALL DEMES WITH SAME SAL TOgether? **WHAT IS THIS**
* `cor_temp` correlation of mutation allele frequency and temperature at end of simulation (all individuals) - IS THIS CALCULATED BY EACH DEME OR ALL DEMES WITH SAME SAL TOgether? *WHAT IS THIS**
* `mutSalEffect` Effect of mutation on salinity
* `mutTempEffect` Effect of mutation on temperature

### 2549986039929_plusneut_MAF01.recode2.vcf.gz

* The VCF file output after pyslim and filtered to an MAF of 0.01
* Each row is a locus simulated in SLiM or added by pyslim
* Each column is an individual (all output)
* The ID of the causal mutations is in the ALT table and should match `mutID` in other tables

### 2549986039929_VCF_causal.vcf.gz

The VCF file output 

* Each row is a locus simulated in SLiM
* Each column is an individual (all output)
* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
* This table is redundant, but used to cross-check the pyslim outputs.

### 2549986039929.trees

The trees file output by SliM that is piped into pyslim.

### 2549986039929_muts_analysis.txt

After each mutation is analyzed in R, 

SEE APRIL 2021 NOTEBOOK POSTS

### 2549986039929_sim_analysis.txt

After each simulation is analyzed in R, 

SEE APRIL 2021 NOTEBOOK POSTS
