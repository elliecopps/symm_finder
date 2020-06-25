"""
energies.py is a vague name but it's an artifact of what this file used to be for and now it's imported by a bunch of
other files so I don't want to change the name

Prints out information about an intance/seed, including degeneracies and swap symmetries plus some other stuff

by default loads a Jij matrix from G_Matrices folder, but can also generate a random Jij from seed
"""

import tfim
import numpy as np
import itertools as it
import tfim
import SwapSymmetry as SS
import argparse

def main():


	parser = argparse.ArgumentParser(description=('Prints degeneracies and swap symmetries'))
	parser.add_argument('-N', type = int, default = 8, help = 'System size')
	parser.add_argument('-p', type = int, default = 0, help = 'instance/seed to print')
	parser.add_argument('-rand_Jij', action = 'store_true', help = "Generate Jij from seed")
	args = parser.parse_args()


	N = args.N
	permutation = args.p
	rand_Jij = args.rand_Jij


	lattice = tfim.Lattice([N])
	basis = tfim.IsingBasis(lattice)

	if rand_Jij:
		G = tfim.Jij_instance(N,1,"bimodal",permutation, True)
	else:
		G = np.loadtxt("G_Matrices/" + str(N) + "/" + str(permutation) + "-G.dat")
		
	# print(G)

	print(G)
	Jij = makeJij(G,N)
	print(Jij)


	ea = energyArray(basis,N,Jij)
	# print("energy array")
	# print(ea)

	print("N: ")
	print(N)
	print("Permutation: " + str(permutation))

	gs = findGroundStates(ea)
	print("min energy:")
	print(ea[gs[0]])
	print("ground state indeces:")
	print(gs)
	print("ground states:")
	print(stateArray(basis,gs))

	pList = pairList(basis,gs)
	print("pair list:")
	print(pList)

	print("overlap: ")
	overlap(basis,gs)

	hCList = hammClusterList(basis,pList)
	print("Hamming cluster list:")
	print(hCList)

	swapsyms = SS.FindSymmetries(N,permutation)
	print("swap symmetries:")
	print(swapsyms)

	swapCList = swapClusterList(basis,swapsyms,pList)
	print("swap symmetry cluster list:")
	print(swapCList)

	aSwapCList = antiSwapClusterList(basis,swapsyms,pList)
	print("\"anti\" swap symmetry cluster list:")
	print(aSwapCList)




	swapclusters = allSwapCLusters(swapCList,aSwapCList)
	print("swap symmetry clusters:")
	print(swapclusters)


def makeJij(G,N):
	"""e
	Turns Jij matrix from form built in tfim.py to standard Jij where J[i][j] is the bond between spins i and j
	"""
	Jij = np.zeros((N,N))
	for j in range(N//2):
		for i in range(N):
			Jij[i][(i-j+N-1) % N] = Jij[(i-j+N-1) % N][i] = G[j][i]
	return Jij

def getEnergy(N,Jij,config):
	Energy = 0
	for i in range(N):
		for j in range(N-i):
			j+=i
			Energy += -Jij[i][j]*(2*config[i]-1)*(2*config[j]-1)
	return Energy

def energyArray(basis,N,Jij):
	Earray = [0 for i in range(basis.M)]
	for i in range(basis.M):
		config = list(map(int,list(bin(i)[2:].zfill(N))))
		Earray[i] = getEnergy(N,Jij,config)

	return Earray

def findSeed(N,Gfind):
	for seed in range(10000):
		G = tfim.Jij_instance(N,1,"bimodal",seed)
		if np.array_equal(G,Gfind):
			print(G)
			print(Gfind)
			print(seed)
			break

def checkSymmetry(N, Jij):
	for i in range(2**N):
		print(i)
		config = list(bin(i)[2:].zfill(N))
		e1 = getEnergy(N,Jij,list(map(int,config)))
		config[2],config[3] = config[3],config[0]
		e2 = getEnergy(N,Jij,list(map(int,config)))

		print(e1,e2)


def findGroundStates(ea):
	mnm = min(ea)

	gl = []
	for i in range(len(ea)):
		if ea[i] == mnm:
			gl.append(i)

	return gl

def minHammDist(basis,s1,s2):
		mDist = np.inf
		for g1i in [i1 for sl1 in s1 for i1 in sl1]:
			for g2i in [i2 for sl2 in s2 for i2 in sl2]:
				dist = 0
				for i in range(len(basis.state(g1i))):
					if basis.state(g1i)[i] != basis.state(g2i)[i]:
						dist += 1
				if dist<mDist:
					mDist = dist
		return mDist


def stateArray(basis,gs):
	sa = []
	for g in gs:
		sa.append(basis.state(g))
	return np.array(sa)

def pairList(basis,groundstates):
	gs = [g for g in groundstates]
	if(len(gs) < 2):
		return [[gs]]
	else:
		pList = []
		while(len(gs) != 0):
			for i in range(1,len(gs)):
				if np.array_equal(np.array(basis.state(gs[0]) + basis.state(gs[i])),np.array([1 for j in range(basis.N)])):
					pList.append([gs[0],gs[i]])
					gs.remove(gs[i])
					gs.remove(gs[0])
					break
		return pList

def hammClusterList(basis,pList):
	hCList = [[p] for p in pList]
	if len(pList) == 1:
		return hCList
	j = 0
	while(True):
		for i in range(j+1, len(hCList)):
			if minHammDist(basis, hCList[j+0], hCList[i]) <= 2:
				hCList[0]+=hCList[i]
				hCList.remove(hCList[i])
				j -= 1
				break

		j+=1
		if j == len(hCList) - 1:
			break

	return hCList 

def swapSym(basis,p1,p2,swaptuple):
	for s1 in p1:
		for s2 in p2:
			flipState = basis.state(s1)
			basis.flip(flipState,swaptuple[0])
			basis.flip(flipState,swaptuple[1])
			if np.array_equal(flipState,basis.state(s2)):
				return True
	return False

def swapClusterList(basis,swapsyms,pList):
	swapCList = {sp:[] for sp in swapsyms[0]}

	for i in range(len(pList)):
		pL = pList[i:]
		for j in range(1,len(pL)):
			for sp in swapsyms[0]:
				if swapSym(basis,pL[0],pL[j],sp):
					swapCList[sp].append([pL[0],pL[j]])
	return swapCList

def antiSwapClusterList(basis,swapsyms,pList):
	aSwapCList = {sp:[] for sp in swapsyms[1]}

	for i in range(len(pList)):
		pL = pList[i:]
		for j in range(1,len(pL)):
			for sp in swapsyms[1]:
				if swapSym(basis,pL[0],pL[j],sp):
					aSwapCList[sp].append([pL[0],pL[j]])
	return aSwapCList


def allSwapCLusters(swapCList,aSwapCList):

	swapclusters = []
	allSwapsList = swapCList.copy()
	allSwapsList.update(aSwapCList)
	for d in allSwapsList:
		if len(allSwapsList[d]) != 0:
			for s in allSwapsList[d]:
				swapclusters.append([tuple(s[0]),tuple(s[1])])

	if len(swapclusters) <= 1:
		return swapclusters



	j = 0
	while(True):
		i = 1
		while(True):
			if i != j:
				if len(set(swapclusters[j]) & set(swapclusters[i])) > 0:
					swapclusters[j] = list(set(swapclusters[j]) | set(swapclusters[i]))
					swapclusters.remove(swapclusters[i])
					i -= 1
			i+=1

			if i == len(swapclusters):
				break

		if j+1 == len(swapclusters):
			break
		j+=1

	return swapclusters


def overlap(basis,gs,prnt=True):

	gspaircombos = list(it.combinations(gs,2))

	minov = np.inf

	for c in gspaircombos:
		overlap = 0
		state0 = basis.state(c[0])
		state1 = basis.state(c[1])
		for i in range(len(state0)):
			if state0[i]==state1[i]:
				overlap+=1

		if overlap < minov and overlap != 0:
			minov = overlap
		if prnt:
			print(str(c) + " " + str(overlap))

	if prnt:
		print(minov)



	return minov



if __name__ == "__main__":
	main()