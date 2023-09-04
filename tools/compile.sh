#!/bin/bash

cd .
python3 setup.py build_ext --inplace
cd py_qcprot
python3 setup.py build_ext --inplace
cd ..