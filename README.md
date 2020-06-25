# symm_finder
The sym_mod file contains all of the functions used to find symmetries in systems of spins.

Running Symmetry_Finder.py from the command line will return symmetries broked into the categories swap, anti_swap, and other
  Arguments: N, seed, # of seeds to test beyond the first one, max_energy
  N = number of spins in system
  seed = the seed corresponding to a particular Jij matrix
  The # of seeds specifies how many consecutive seeds you want to test for symmetries
  A max_energy of 0 corresponds to checking all energies, 1 corresponds to checking only the ground state

sym_mod.pyc is a Jupyter Notebook with all of the functions that I use.
