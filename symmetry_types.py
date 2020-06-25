import tfim
import numpy as np
import itertools as it
import energies as en
import argparse
import easydict



'''
Code determines whether a seed or group of seeds has swap or anti-swap symmetries.
Modified from Asher's symcheck.py
'''


def main(seed, N):
    
    args = easydict.EasyDict({'cs':2})

    clustersize = args.cs

    seeds = [seed]

    for seed in seeds:
        symarray = []
        swap_symmetries = []
        anti_swap_symmetries = []

        symsy = 0
        for clustersize in range(1, N-2):
            if clustersize == N/2:
                continue

            G = tfim.Jij_instance(N,1,"bimodal",seed,True)

            Jij = makeJij(G,N)


            combinations= list(it.combinations([i for i in range(N)], clustersize*2))
            
            for c in combinations:
				#this code below is to make sure all combinations of groups of clusters are checked in each
				#possible combination of spins
                
                paircombos = list(it.combinations(c,clustersize))
                pairlist = [paircombos[i] + paircombos[len(paircombos) - i -1] for i in range(len(paircombos)//2)]

                for p in pairlist:

                    
                    remove_index = list(p)
                    
                    anti_swap = True
                    swap = True
                    
                    for ind in range(len(p)//2):
                        vprime = Jij[p[2*ind]].copy()
                        compare = Jij[p[2*ind + 1]].copy()
                        for ind2 in range(len(p)//2):
                            vprime[p[2*ind2]], vprime[p[2*ind2 + 1]] = vprime[p[2*ind2 + 1]], vprime[p[2*ind2]]
                        
                        vprime = np.delete(vprime, remove_index)
                        
                        compare = np.delete(compare, remove_index)
                        
                        
                        if (np.array_equal(compare,vprime)):
                            anti_swap = False
                            continue
                        elif (np.array_equal(compare, -vprime)):
                            swap = False
                            continue
                        else:
                            swap = False
                            anti_swap = False
                            
                    if swap:
                        symsy = clustersize

                        swap_symmetries.append(p)
                    
                    if anti_swap:
                        symsy = clustersize
                        anti_swap_symmetries.append(p)

    return swap_symmetries, anti_swap_symmetries



def makeJij(G,N):
	"""
	Turns Jij matrix from form built in tfim.py to standard Jij where J[i][j] is the bond between spins i and j
	"""
	Jij = np.zeros((N,N))
	for j in range(N//2):
		for i in range(N):
			Jij[i][(i-j+N-1) % N] = Jij[(i-j+N-1) % N][i] = G[j][i]
	return Jij



