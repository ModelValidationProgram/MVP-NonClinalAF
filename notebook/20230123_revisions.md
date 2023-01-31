## 20230120

A number of tweaks were made to the tutorial and the multivariate
simulation code, which are documented in the response to reviewers.

Cross-checking all the metadata and code is synced to github:

- [x] sync github with laptop - make sure no lingering changes on the HPC

2023-01-26 
```
git status

#       new file:   multipheno_multienvi/SlimMultvar20221024.err
#       new file:   multipheno_multienvi/SlimMultvar20221024.out
#       modified:   multipheno_multienvi/output_multisim/892657863_ind.txt
#       modified:   multipheno_multienvi/output_multisim/_pytree.out.txt
#       modified:   src/g-CheckRunsAreComplete.Rmd
#       modified:   src/g-FinalAnalysis.Rmd
#       new file:   src/g-FinalAnalysis.html
#       new file:   src/g-FinalAnalysis.pdf
#       modified:   summary_20220428_20220726.txt
#       modified:   summary_20220428_20220726_RDApredictions.txt
```
I synced everything today and cleaned up the working directories. 


- [x] Tutorial
	- [x] check that someone can reproduce it
- [x] `README.md` - revisit this after deciding whether to keep the figures and example folder
- [x] `summary_20220428_20220726_RDApredictions.txt_metadata.md`
- [x] `summary_20220428_20220726.txt_metadata.md`
		- checked the metadata and the outputs
- [x] `summary_LAcalc_20220428_20220726.txt_metadata.md`
- [x] `summary_LAthroughtime_20220428_20220726.txt_metadata.md`

- [x] `/sim_ouput_20220428_metadata/` 
- [x] README.md
- [x] seed_ind.txt_metadata.md
- [x] seed_info.txt_metadata.md
- [x] seed_LA.txt_metadata.md
- [x] seed_muts.txt_metadata.md
- [x] seed_popInfo_m.txt_metadata.md
- [x] seed_popInfo.txt_metadata.md
- [x] seed_Rout_af_pop.txt_metadata.md
- [x] seed_Rout_af_sal.txt_metadata.md
- [x] seed_Rout_af_temp.txt_metadata.md
- [x] seed_Rout_ind_subset.txt_metadata.md
		- I checked it, but double check the numbering on $opt_0$ in the newer outputs
- [x] seed_Rout_af_pop.txt_metadata.md
- [x] seed_Rout_muts_full.txt_metadata.md
- [x] seed_Rout_RDA_predictions.txt_metadata.md


 
 
- [x] `multipheno_multienvi/`

- [x] README

- [x] multipheno_multienvi/bioclim
	- [x] adaptive_env.txt
	- [x] bioclim.txt
	- [x] MAT_BC_360x360.csv
	- [x] MTDQ_BC_360x360.csv
	- [x] MTWetQ_BC_360x360.csv
	- [x] nuisance_env.txt
	- [x] PDM_BC_360x360.csv
	- [x] PWarmQ_BC_360x360.csv
	- [x] PWM_BC_360x360.csv
	- [x] README

- [x] multipheno_multienvi/output_multisim
	- [x] README.md metadata for all files in the folder here

- [x] multipheno_multienvi/src
	- [x] a-BIOCLIM_extraction.R
	- [x] b-non_wf_range_exp_working.slim
	- [x] c-pyslim.py
	- [x] c-pyslim.sh
	- [x] d-processVCF.Rmd --> RERUN ON LAPTOP AND MAKE PDF
	- [x] README.md


- [x] re-run the multivariate code to make sure my cleaning up didn't break anything



- [ ] `src/` 
	- [ ] README.md
	- [ ] env
	- [ ] 0a-setUpSims.Rmd
	- [ ] 0b-final_params-20220428.RData
	- [ ] 0b-final_params-20220428.txt
	- [ ] 0b-final_params-20220428.txt_metadata.md
	- [ ] 0b-final_params-fastruns-20220428-b.txt
	- [ ] 0b-final_params-fastruns-20220428.txt
	- [ ] 0b-final_params-longruns-20220428.txt
	- [ ] a-PleiotropyDemog_20220428.slim
	- [ ] b-process_trees.py
	- [ ] c-AnalyzeSimOutput.R
	- [ ] c2-ContributionToLA.R
	- [ ] d-run_nonAF_sims_0Slim-fastruns-20220428-b.sh
	- [ ] d-run_nonAF_sims_0Slim-fastruns-20220428.sh
	- [ ] d-run_nonAF_sims_0Slim-longruns-20220428.sh
	- [ ] e-run_nonAF_sims_1R-fastruns-20220428-b.sh
	- [ ] e-run_nonAF_sims_1R-fastruns-20220428.sh
	- [ ] e-run_nonAF_sims_1R-longruns-20220428.sh
	- [ ] f-catoutputs.sh
	- [ ] g-CheckRunsAreComplete.Rmd
	- [ ] g-FinalAnalysis.html
	- [ ] g-FinalAnalysis.pdf -- DOWNLOAD AND KEEP A COPY TO ENSURE REPROD.
	- [ ] g-FinalAnalysis.Rmd -- CHECK TO MAKE SURE RUNS ON CLUSTER AND MAKE PDF
	- [ ] g2-FinalAnalysis-ContributionToLA.Rmd
	- [ ] g3-LAthroughTime.Rmd
	- [ ] g4-DifferentArchSameGenotype.Rmd
	- [ ] parameterOrderForShellScript.txt
	- [ ] plot_mig_Graphs.R

