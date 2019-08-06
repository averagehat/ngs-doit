#!/bin/bash
set -e
set -x

if [[ -z "$1" ]]; then
    echo "Requires install prefix for new miniconda as an argument. i.e. /home/user/vartable/miniconda"
fi 
    
curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda3.sh
bash miniconda3.sh -p "$1" -b 

export PATH=$1/bin/:$PATH

if [ "$PYTHON_VERSION" != "3.7" ]
    then conda install -y python="$PYTHON_VERSION"
fi

conda update -y -n base -c defaults conda # want conda 4.7

conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y --file requirements-conda.txt

pip install -r pip-requirements.txt
pip install git+https://github.com/python/mypy   # currently requires mypy for typing_extensions
python setup.py install
