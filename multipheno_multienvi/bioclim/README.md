- `multipheno_multienvi/bioclim`

This Directory contains the bioclim environments used as selective environments
in the continuous space simulation

- `adaptive_env.txt` This file describes the environmental value for selective environments 
at each long,lat location extracted from the BIOCLIM database
	- `x` longitude
	- `y` latitude
	- `MAT` environmental value for this variable
	- `MTWetQ` environmental value for this variable
	- `MTDQ` environmental value for this variable
	- `PWM` environmental value for this variable
	- `PDM` environmental value for this variable
	- `PWarmQ` environmental value for this variable
	- `slim_x` x location in slim
	- `slim_y` y location in slim
	
- `nuisance_env.txt` This file describes the environmental value for non-selective environments
at each long,lat location extracted from the BIOCLIM database
	- `x` longitude
	- `y` latitude
	- `ISO` environmental value for this variable
	- `TSsd` environmental value for this variable
	- `PSsd`environmental value for this variable
	- `slim_x` x location in slim
	- `slim_y` y location in slim
	
- `bioclim.txt`
	- `BIO` BIOCLIM variable
	- `ABBRV` Abbreviation used in this study
	- `DESC` Description of the variable

The following 6 files were used as input to the `slim` code. 
Each file describes the environmental value at that location on a 360 x 360 grid.

- `MAT_BC_360x360.csv`
- `MTDQ_BC_360x360.csv`
- `MTWetQ_BC_360x360.csv`
- `PDM_BC_360x360.csv`
- `PWarmQ_BC_360x360.csv`
- `PWM_BC_360x360.csv`