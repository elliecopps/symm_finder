#!/usr/bin/env python
# coding: utf-8


import sym_mod
import symmetry_types
import argparse
import tfim

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int, help='Number of spins')
parser.add_argument('seed', type=int, help='Jij matrix seed')
parser.add_argument('max_energy', type=int, help='Maximum energy level to search for symmetries in')
args = parser.parse_args()
    
N = args.N
seed = args.seed
max_energy = args.max_energy


PBC = True
L = [N]
lattice = tfim.Lattice(L, PBC)
J = 1

Jij = tfim.Jij_instance(N,J,dist="bimodal",seed=seed,even=True)
basis = tfim.IsingBasis(lattice)


grouped_configurations = sym_mod.grouped_configs(basis, Jij)


transformations = sym_mod.list_transformations(grouped_configurations, max_energy, N, basis)


symmetries = sym_mod.symmetry(transformations, grouped_configurations, basis)


sorted_list = sym_mod.list_sorter(symmetries, N)



simplified_sym = sym_mod.remove_combo(sorted_list)



sorted_sym = symmetry_types.main(seed, N)


classified_sym = sym_mod.symmetry_sorter(simplified_sym, sorted_sym)

#if __name__ == "__main__":
    #print classified_sym
