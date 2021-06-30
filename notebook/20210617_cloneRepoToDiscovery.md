# 20210617 

Cloned the github repo to Discovery

First, I followed these instructions to set up a public/private key connection between GitHub and Discovery: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh 

```
ssh lotterhos@discovery.neu.edu
# signed in
cd work/lotterhos
git clone git@github.com:ModelValidationProgram/MVP-NonClinalAF.git

git branch alan
git checkout alan
git commit -a -m "add branch"
```

## To do for Alan:
- [x] make sure can access work/lotterhos/MVP-NonClinalAF
- [x] make sure can link github account to discovery and make changes on `alan` branch
- [x] don't "F" up using github - review the "git" chapter in Buffalo's book and code below
- [x] review the command line script in `run_nonAF_sims_alan.sh`. 
- [x] Create a conda environment for running the shell script
    - [x] check slim is verson 3.6
    - [x] the python line takes 24 hours on my laptop
    - [x] update the notebook and push to github
- [x] test the bash script gives proper output, in particular the VCF editing line
    - [x] update the notebook and push to github
- [ ] Create a bash script for running 1 replicate of the simulation as shown in the bash script
    - [ ] update the notebook and push to github
- [ ] review the code in `0c-writeSimScript.R` that creates some of the scripts
- [ ] Write an array to run the first 150 lines (parameter combos) of the `0b-final_params.txt`
    - [ ] update the notebook and push script to github. Note here there are some simulation outputs that are not pushed to github because they are so big. So we should review those before you try to push the outputs too.

#### Step 1: get off the login node

srun

#### Step 2 view the branches:

`git branch` 

alan

* main

#### Swith to branch `alan`:

`git checkout alan` 

Switched to branch 'alan'


#### Tracking files

Let’s use git add to tell Git to track the files:
* `git add README data/README` will track both README in the current dir and data/README
* `git add --all sim_outputs1/` adds everything in the folder

`git status` tells us what is staged to be committed

after adding everything, we commit it `git commit -a -m "reogoranize"`

then you can push `git push origin alan` for alan's branch, or  `git push origin main` for the main branch

Once we are happy with changes on alan's branch, we will submit a merge request with the main branch
