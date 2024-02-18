GUI screen:
1. open the terminal:
ssh ram -J huanliu@raptor.acm.unc.edu -L 5555:localhost:3389
2. open the RDP in the laptop and connect to 'localhost:5555'

log in server:
1. connect the server
2. ssh ram -J huanliu@raptor.acm.unc.edu
(PIN:dgybis888666#)
3. cd /ram/USERS/huanliu/1_filtered_wavelet/

the guidence of using gspbox:
1. start matlab: matlab (needed activated)
2. gsp_start
3. gsp_install
4. add ./utils to Path

project code order:
I. produce the wavelet:
Filter_identify.m
II. main procedure:
1. Generate_CFC(0,2,7);
2. Linear_Regression(2,7,1);
III. post process:
1. Mat_to_Txt(2,7,1);
2. HCP_MMP_binary_GLM_on_CFC(2,7,3,360);

conda: cd /export_home/huanliu/
. ./activate
conda config --set auto_activate_base false
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create -n R-3.5.1
conda info --envs
source activate R-3.5.1
conda install r-base=3.5.1
conda remove --name R-3.5.1 --all
conda deactivate

Rscript LR.R
