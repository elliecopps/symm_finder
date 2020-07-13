# symm_finder

The sym_mod.py file contains all of the functions that are used in the main Symmetry_Finder.py file

From the command line, python Symmetry_Finder.py will run the function. The command line arguments are N, seed, # of seeds to test, and maximum energy level. To test only the ground state, max_energy should be 1. 0 will test all energy levels.

Currently, the function only returns symmetries and transformations with a hamming distance less than N/2. It also throws out any transformations of hamming distance 4 that consist of two transformations of hamming distance 2. In the future, we may add functionality to throw out higher-order combinations and to find symmetries with greater hamming distances.

Symmetry_Finder.py writes data to three files. degeneracies.json contains a list of lists with seeds and clusters of ground states that are related by symmetries. symmetries.json has a list of seeds and all of the symmetries associated with them grouped into swap, anti swap, and other.

Must be running python 2.
