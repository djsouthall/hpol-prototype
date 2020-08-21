# rnog-hpol-prototype
This contains code and data created while working towards a finalized hpol antenna for RNO-G.

I use some functions that are defined from another package: [plot_ff.py](https://github.com/djsouthall/beacon/blob/master/tools/plot_ff.py).  These aren't complicated and could easily be redifined in this project, but I am the only one using it and this is fine for me. 

To use certain tools and scripts in this repository you may need to define the environment variables *BEACON_ANALYSIS_DIR* (as specified in the BEACON analysis code framework [BEACON](https://github.com/djsouthall/beacon)), as well as *HPOL_ANALYSIS_DIR*, which is similar but pointing to the folder containing this repository on your local machine. 

Additionally I am moving towards using direct imports of packages, which will require you to add these paths to your python path.

For example in my bashrc I have:
HPOL_ANALYSIS_DIR="/home/dsouthall/Projects/Greenland/hpol_prototype/"
export HPOL_ANALYSIS_DIR
export PYTHONPATH=$PYTHONPATH:/home/dsouthall/Projects/Greenland/
