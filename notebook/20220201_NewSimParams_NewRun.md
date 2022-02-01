
# Revise Simulation parameters to give more realistic FST

### Edit `0a-setUpSims.Rmd`
Old simulation parameters:
```
c("SS-Clines", "SS-Mtn", "Est-Clines")
MIG_x = c(0.01, 0.01, 0.49)
MIG_y = c(0.01, 0.01, 0.01)
```

New simulation parameters in `0a-setUpSims.Rmd`:
```
MIG_x = c(0.1, 0.1, 0.49)
MIG_y = c(0.1, 0.1, 0.05)
```

(Need to run `0a-setUpSims.Rmd` with SimID 20220201)

### New simulation parameters in `a-PleiotropyDemog_20220101.slim`:

The m-breaks scenario and variable-m scenario:

Old code:
```
	if (demog=="Estuary"){

	var_m_estuary = c(0.001,rep(0.01,10),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.001;
```

New code:
```
	var_m_estuary = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.01;
```

Old code:
```
	if (demog=="SS"){ // STEPPING STONE DEMOGRAPHY

	var_m_ss = c(rep(0.001,6),0.01,0.1,0.25);
  # mig_breaks was 0.0001
```  

New code:
```
  var_m_SS = c(0.005, rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_SS = 0.01;
```

### Updating pipeline files

Updated files for a 20220201 simulation ID

```
(base) [lotterhos@login-01 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201.sh
Submitted batch job 23002489
```
