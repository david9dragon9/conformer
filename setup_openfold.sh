#!/bin/bash
git clone https://github.com/aqlaboratory/openfold.git
rm openfold/openfold/np/relax/*
touch openfold/openfold/np/relax/relax.py
touch openfold/openfold/np/relax/__init__.py
cd openfold
python3 setup.py install
cd ..
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
