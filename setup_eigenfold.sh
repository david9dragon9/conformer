#!/bin/bash
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.0.0+cu118.html
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cu118.html
git clone https://github.com/bjing2016/EigenFold.git
wget https://helixon.s3.amazonaws.com/release1.pt
git clone https://github.com/bjing2016/OmegaFold.git OmegaFoldE