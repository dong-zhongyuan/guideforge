# RibonanzaNet Inference Wrapper
Self-contained, automatable inference for the RibonanzaNet series of models

## Setup
* Install conda
* Clone this repository (make sure to use `--recursive`, or after cloning run `git submodule update --init`)
* `conda env create -f RibonanzaNet/env.yml -p .venv`
* `.venv/bin/pip install -r requirements.txt`

## Usage
```sh
.venv/bin/python src/rnet-shape.py GGGGAAAACCCC
.venv/bin/python src/rnet-2d.py GGGGAAAACCCC
.venv/bin/python src/rnet-deg.py GGGGAAAACCCC
.venv/bin/python src/rnet-drop.py GGGGAAAACCCC
```
