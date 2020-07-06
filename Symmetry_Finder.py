#!/usr/bin/env python
# coding: utf-8


import sym_mod
import symmetry_types
import argparse
import tfim
import cProfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('N', type=int, help='Number of spins')
    parser.add_argument('seed', type=int, help='Jij matrix seed')
    parser.add_argument('seed_range', type=int, help='Number of seeds to test')
    parser.add_argument('max_energy', type=int, help='Maximum energy level to search for symmetries in')
    args = parser.parse_args()
        
    N = args.N
    seed = args.seed
    seed_range = args.seed_range
    max_energy = args.max_energy

    PBC = True
    L = [N]
    lattice = tfim.Lattice(L, PBC)
    J = 1
    basis = tfim.IsingBasis(lattice)

    for seed in range(seed, (seed+seed_range)):

        print 'Seed: ', seed

        Jij = tfim.Jij_instance(N,J,dist="bimodal",seed=seed,even=True)
        
        Energies = -1 * (tfim.JZZ_SK_ME(basis,Jij))

        grouped_configurations = sym_mod.grouped_configs(Energies)

        transformations = sym_mod.list_transformations(grouped_configurations, max_energy, N, basis)

        symmetries = sym_mod.symmetry(transformations, grouped_configurations, basis)
        
        if len(symmetries) == 0:
            continue

        sorted_list = sym_mod.list_sorter(symmetries, N)

        simplified_sym = sym_mod.remove_combo(sorted_list, basis)

        sorted_sym = symmetry_types.main(seed, N)
        #print sorted_sym
        
        classified_sym = sym_mod.symmetry_sorter(simplified_sym, sorted_sym)
        
        print ' '
        

#cProfile.run('main()')
if __name__ == '__main__':
    main()
