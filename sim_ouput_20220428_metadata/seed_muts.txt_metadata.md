
### seed_muts.txt (SLiM output) mutation table
This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99

* `seed`  simulation seed
* `mutID` mutation ID - should match mutID in other tables and VCF file
* `muttype` mutation type in SLiM
* `p` allele frequency of derived mutation - not actually the minor allele frequency
* `cor_sal` correlation of mutation allele frequency and salinity at end of simulation (all individuals) - each deme is a datapoint
* `cor_temp` correlation of mutation allele frequency and temperature at end of simulation (all individuals) - each deme is a datapoint
* `mutSalEffect` Effect of mutation on salinity
* `mutTempEffect` Effect of mutation on temperature