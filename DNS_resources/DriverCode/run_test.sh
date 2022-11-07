#!/bin/bash

# Author: Radu Cimpeanu
# Date: 04/11/2022

# Copy the MasterImpact folder to a test folder
cp -r MasterImpact/ Impact_Test
cd Impact_Test/		

# Compile the code (video capabilities enabled)
qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11

# Tailor the CPU count to desired resources
export OMP_NUM_THREADS=4

# The parameters passed on to the code are:
# 1: Weber number
# 2: Dimensionless pool height
# 3: Maximum resolution level
# 4: Maximum runtime (dimensionless)
# Details on the definitions above can be found in the manuscript
	
# Execute the code
./DropImpact 345 0.3 9 0.501
