# RNO-G HPol Antenna Design

This design was finalized (for the current deployment cycle), and several dozen antennas have been built with this design.  Images of the finalized design can be seen below.  The antennas consist of 8in OD alluminum tubes with custom machining for slots and structural support holes.  The centerfeed consists of a connected custom PCB circuits, supported with several waterjet nylon inserts.  Each end of the antennas have waterjet cut nylon caps, with the top having additional stucture to hold the front-end electronics.  These endcaps have rope strung through which is ultimately used to hang the antenna in the bore-hole.  

# Images
Opaque Element             |  Transparent Element
:-------------------------:|:-------------------------:
<img src="https://github.com/djsouthall/hpol-prototype/blob/master/Dsouthall-Hpol-figure_opaque_4k.png?raw=true" alt="Opaque Element" width="500"/>  |  <img src="https://github.com/djsouthall/hpol-prototype/blob/master/Dsouthall-Hpol-figure_transparent_4k.png?raw=true" alt="Transparent Element" width="500"/>

Feed
:--------------------------------------------------:
| <img src="https://github.com/djsouthall/hpol-prototype/blob/master/Dsouthall-assembled-feed_4k.png?raw=true" alt="Antenna Feed" width="1000"/> |



# Development Notes
This contains code, data, and design files created while working towards a finalized hpol antenna for RNO-G.

I use some functions that are defined from another package: [plot_ff.py](https://github.com/djsouthall/beacon/blob/master/tools/plot_ff.py).  These aren't complicated and could easily be redifined in this project, but I am the only one using it and this is fine for me. 

To use certain tools and scripts in this repository you may need to define the environment variables *BEACON_ANALYSIS_DIR* (as specified in the BEACON analysis code framework [BEACON](https://github.com/djsouthall/beacon)), as well as *HPOL_ANALYSIS_DIR*, which is similar but pointing to the folder containing this repository on your local machine. 

Additionally I am moving towards using direct imports of packages, which will require you to add these paths to your python path.

For example in my bashrc I have:
HPOL_ANALYSIS_DIR="/home/dsouthall/Projects/Greenland/hpol_prototype/"
export HPOL_ANALYSIS_DIR
export PYTHONPATH=$PYTHONPATH:/home/dsouthall/Projects/Greenland/
