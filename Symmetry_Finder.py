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
    
    f = open("symmetries.txt", "w+")
    g = open("degeneracies.txt", "w+")

    for seed in range(seed, (seed+seed_range)):


        print 'Seed: ', seed
        
        f.write('Seed: ' + str(seed) + '\n')
        g.write('Seed: ' + str(seed) + '\n')

        Jij = tfim.Jij_instance(N,J,dist="bimodal",seed=seed,even=True)
        
        Energies = -1 * (tfim.JZZ_SK_ME(basis,Jij))

        grouped_configurations = sym_mod.grouped_configs(Energies)
        
        g.write("Ground state degeneracy: " + str(grouped_configurations[0]) + '\n')
        
        transformations = sym_mod.list_transformations(grouped_configurations, max_energy, N, basis)

        symmetries = sym_mod.symmetry(transformations[0], grouped_configurations, basis)
        
        if len(symmetries) == 0:
            f.write('No symmetries' + '\n')
            f.write('\n')
            g.write('Degeneracy pairs related by symmetries: ' + '\n')
            g.write('     ' + 'None' + '\n')
            g.write('\n')
            
            continue

        sorted_list = sym_mod.list_sorter(symmetries, N)

        simplified_sym = sym_mod.remove_combo(sorted_list, basis)

        sorted_sym = symmetry_types.main(seed, N)
        
        classified_sym = sym_mod.symmetry_sorter(simplified_sym, sorted_sym)
        
        if max_energy == 1:
            degen_clusters = sym_mod.degen_sym(transformations[1], classified_sym)
            
            g.write('Degeneracy pairs related by symmetries: ')
            if len(degen_clusters) == 0:
                g.write('None' + '\n')
            else:
                for cluster in degen_clusters:
                    g.write('\n' + '     ' + str(cluster[0]) + ' , ' + str(cluster[1]) + '\n')
                    g.write('     ' + 'Transformation: ' + str(cluster[2]) + '\n')
                    g.write('\n')
        
        f.write('Swap symmetries: ' + str(classified_sym[0]) + '\n')
        f.write('Anti-swap symmetries: ' + str(classified_sym[1]) + '\n')
        f.write('Other: ' + str(classified_sym[2]) + '\n')
        f.write('\n')
        
        print ' '
        
    f.close()
    g.close()

#cProfile.run('main()')
if __name__ == '__main__':
    main()
