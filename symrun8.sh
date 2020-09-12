#!/bin/bash

k=$SLURM_ARRAY_TASK_ID

python Symmetry_Finder.py -N 8 -seed 0 -seed_range 1000 -max_energy 1