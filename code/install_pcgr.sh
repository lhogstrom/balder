PCGR_VERSION="1.1.0"
PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
PLATFORM="linux" 

# mamba is a much faster alternative to conda
/home/larsonh/miniconda3/condabin/conda install mamba -c conda-forge

/home/larsonh/miniconda3/condabin/mamba create --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock --prefix ./pcgr
/home/larsonh/miniconda3/condabin/mamba create --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock --prefix ./pcgrr

# you need to specify the directory of the conda env when using --prefix
/home/larsonh/miniconda3/condabin/:qconda activate ./pcgr

# test that it works
pcgr --version
