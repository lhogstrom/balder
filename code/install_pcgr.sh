PCGR_VERSION="1.2.0"
PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
PLATFORM="linux" 

cd /data/larsonh/other

# mamba is a much faster alternative to conda
/home/larsonh/miniconda3/condabin/conda install mamba -c conda-forge

/home/larsonh/miniconda3/condabin/mamba create --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock --prefix ./pcgr
/home/larsonh/miniconda3/condabin/mamba create --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock --prefix ./pcgrr

source /home/larsonh/miniconda3/etc/profile.d/conda.sh

# you need to specify the directory of the conda env when using --prefix
conda activate ./pcgr

# test that it works
pcgr --version
