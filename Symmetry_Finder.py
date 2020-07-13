#!/usr/bin/env python
# coding: utf-8


import sym_mod
import symmetry_types
import argparse
import tfim
import cProfile
import json

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
    
    f = open("symmetries.txt", "w+") #file containing seeds with the symmetries sorted by type
    g = open("degeneracies.txt", "w+") #contains seeds and ground state clusters
    h = open("groun_sym.txt", "w+") #contains seed and all symmetries of ground states sorted by type
    
    f_list = []
    g_list = []
    h_list = []

    for seed in range(seed, (seed+seed_range)):


        print 'Seed: ', seed

        Jij = tfim.Jij_instance(N,J,dist="bimodal",seed=seed,even=True)
        
        Energies = -1 * (tfim.JZZ_SK_ME(basis,Jij))

        grouped_configurations = sym_mod.grouped_configs(Energies)
        
        transformations = sym_mod.list_transformations(grouped_configurations, max_energy, N, basis)

        symmetries = sym_mod.symmetry(transformations[0], grouped_configurations, basis)
        
        if len(symmetries) == 0:
            continue

        sorted_list = sym_mod.list_sorter(symmetries, N)

        simplified_sym = sym_mod.remove_combo(sorted_list, basis)

        sorted_sym = symmetry_types.main(seed, N)
        
        classified_sym = sym_mod.symmetry_sorter(simplified_sym, sorted_sym)
        symm_data = [seed, classified_sym]
        f_list.append(symm_data)
        
        degen = (sym_mod.degen_sym(transformations[1], grouped_configurations[0], classified_sym))
        vector = sym_mod.sym_vector(degen[0])
        clusters = sym_mod.cluster(vector, grouped_configurations[0])
        
        ground_symm = [seed, degen[1]]
        h_list.append(ground_symm)
        
        seed_w_clusters = [seed,clusters]
        g_list.append(seed_w_clusters)
    
    g.write(json.dumps(g_list))
    f.write(json.dumps(f_list))
    h.write(json.dumps(h_list))
    f.close()
    g.close()
    h.close()

#cProfile.run('main()')
if __name__ == '__main__':
    main()
    
'''
q = open("degeneracies.txt")
data = json.load(q)
print data
'''
#data is a list of lists containing the seed and the ground state clusters


