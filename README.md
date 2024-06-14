# dark_showers_tool

Tool to generate pythia 8 configuration cards for dark shower models. See arXiv 2103.01238 and xxxx.xxxxx for a description of the models and assumptions.

## Requirements
- python 3.6 or higher
- scipy
- numpy

## Included
- dark_shower.py : python package corresponding to the one flavor model in arXiv 2103.01238
- dark_shower_two_flavor.py : python package corresponding to the two flavor model in arXiv xxxx.xxxxx
- data: folder containing lifetime grids which are interpolated by dark_shower.py
- dark_showers_demo.ipynb: example jupyter notebook, which serves as a short tutorial for dark_shower.py
- dark_showers_two_flavor_demo.ipynb: example jupyter notebook, which serves as a short tutorial for dark_shower_two_flavor.py
- higgs_portal.ipynb: example jupyter notebook, which serves as a short tutorial for dark_shower.py
- vector_portal.ipynb: example jupyter notebook, which serves as a short tutorial for dark_shower.py
- card: folder containing a couple example pythia 8 configuration cards, as generated with dark_shower_two_flavor.py

## functionality
- calculates lifetimes and branching ratios for several well motivated models, as described in arXiv 2103.01238 and arXiv xxxx.xxxxx
- provides an estimate of the minimal permissible lifetime for the various decay portals in arXiv 2103.01238
- writes out theoretically consistent pythia 8 configuration files for the models

## additions by Karri
- produce_cards.py # produces a grid of pythia cards
- run_pythia.py # runs pythia


## Single card runner
- make
- ./bin/card_runner.exe pythiaCard outFileName maxEvents

## Run all cards
- The script run_cards.py takes all the cards of the ./cards directory and runs them, creating the hepmc files in the ./output folder. 
- If you want to ignore a card so the run_cards.py doesn't run it, you can add it to the ./cards/IGNORE file