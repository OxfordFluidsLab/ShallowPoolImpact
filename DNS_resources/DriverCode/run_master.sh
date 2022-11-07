#!/bin/bash

# Author: Radu Cimpeanu
# Date: 04/11/2022

# Parameter sweeps
for LEVEL in 11 12 13 14; do
	for WEBER in 345; do
		for HEIGHT in 0.1 0.2 0.3; do

			# Copy the MasterImpact folder to a suitably named clone
			cp -r MasterImpact/ Impact-H$HEIGHT-We$WEBER-Level$LEVEL
			cd Impact-H$HEIGHT-We$WEBER-Level$LEVEL/

			# Compile the code (video capabilities enabled)
			qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
			
			# Tailor the CPU count to desired resources
			export OMP_NUM_THREADS=8

			# The parameters passed on to the code are:
			# 1: Weber number
			# 2: Dimensionless pool height
			# 3: Maximum resolution level
			# 4: Maximum runtime (dimensionless)
			# Details on the definitions above can be found in the manuscript
			
			# Execute the code
			./DropImpact $WEBER $HEIGHT $LEVEL 1.001

			cd ..
		done
	done
done
