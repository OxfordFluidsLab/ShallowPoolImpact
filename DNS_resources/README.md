# Shallow pool impact
## Direct numerical simulation summary
Direct numerical simulation code infrastructure for high-speed drop impact onto liquid pools of varying heights, supporting collaborative work with the [Oxford Fluids Laboratory](https://github.com/OxfordFluidsLab).  

## Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. 
* Full visualisation capabilities have been used in order to generate animations. These may be switched off depending on the local architecture.

## Running the code
Once the Basilisk structure is in place, the driver code here is built in order to navigate parameter sweeps in resolution level, Weber number and dimensionless pool height, with one or several of each values added to the run_master.sh for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. There is also a base version of the running script, called run_test.sh, which is built for debugging purposes and handles a single case.

The code can be executed by simply executing this shell script via *sh run_master.sh* inside a terminal. Output will then be produced within a foldering structure that consists of summary DNS execution information, mass conservation and VOF data, interface coordinates, simulation slices and animations, which can be used for further post-processing.

## Example results
The uploaded framework provides a subset of the data generated for two typical cases, with $We=345$ and two different film heights, $h^* = 0.10$ and $h^* = 0.32$, respectively. These are production runs set to resolution level $13$, which was found to be suitable in terms of generating mesh-independent results for the metrics of interest in our region of interest in the parameter space.
