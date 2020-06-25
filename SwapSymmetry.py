import tfim as tfim
import numpy as np
import itertools as it
import argparse

'''
I don't really use this file anymore but it's imported by energies to find cluster size 1 symmetries
so I'm keeping it here for now.
'''


def makeJij(G,N):
	"""
	Turns Jij matrix from form built in tfim.py to standard Jij where J[i][j] is the bond between spins i and j
	"""
	Jij = np.zeros((N,N))
	for j in range(N//2):
		for i in range(N):
			Jij[i][(i-j+N-1) % N] = Jij[(i-j+N-1) % N][i] = G[j][i]
	return Jij

def FindSymmetries(N,seed):

	lattice = tfim.Lattice([N])
	basis = tfim.IsingBasis(lattice)
	# G = np.loadtxt("G_Matrices/" + str(N) + "/" + str(permutation) + "-G.dat")
	G = tfim.Jij_instance(N,1,"bimodal",seed,True)

	Jij = makeJij(G,N)

	c = list(it.combinations([i for i in range(N)],2))
	tl = []
	tla = []
	for i in range(len(c)):
		s1 = c[i][0]
		s2 = c[i][1]
		if sum(abs(Jij[s1] - Jij[s2])) == 2:
			tl.append((s1,s2))
		if sum(abs(Jij[s1] + Jij[s2])) == 2:
			tla.append((s1,s2))

	return tl, tla



def main():

	parser = argparse.ArgumentParser(description=("Testing for and playing around with swap symmetry in infinite dimensional spin glass") )
	parser.add_argument('-N', type = int, default = 4, help = "Number of spins")
	parser.add_argument('-seed', type = int, default = 1, help = "Seed the randomly generated Jij matrix")
	parser.add_argument('-all', type = bool, default = False, help = "Find all swap symmetries for  every j matrix")

	args = parser.parse_args()
	all = args.all
	N = args.N


	if all:
		for s in range(1,101):
			J_Directory = "J_Matrices/J"+str(N)+"/J_Matrix"+str(s)+"/SwapSymmetries.dat"
			np.savetxt(J_Directory,FindSymmetries(N,s),fmt='%i')



	seed = args.seed
	print(FindSymmetries(N,seed))

	# findSymBps(N,FindSymmetries(N,seed))






	# spectrum = tfim.JZZ_SK_ME(basis,G)
	# energies = np.unique(spectrum)

	# spectrum_array = np.zeros(len(energies))
	# for e in energies:
	# 	print(e)
	# 	for i in range(len(spectrum)):
	# 		if spectrum[i] == e:
	# 			s = basis.state(i)
	# 			print(str(basis.state(i)))


	



def findSymBps(N,tl):	
	pairlist = []
	for t in tl:
		s1 = t[0]
		s2 = t[1]

		combos,permlist = uniqueBps(N)
		symlist = []
		for i in range(len(combos)/2):
			if ((s1 in combos[i]) != (s2 in combos[i])):
				symlist.append(i)
		print(permlist)
		symlistE = list(tuple(symlist))
		while(len(symlist) > 0):
			bp1 = list(permlist[symlist[0]])
			bp2 = standardForm(swapSpins(bp1,s1,s2),permlist)
			ipair = permlist.index(tuple(bp2))
			pair = [symlist[0],ipair]
			pairlist.append(pair)
			symlist.remove(symlist[0])
			symlist.remove(ipair)
	print(pairlist)
	# equivalenceList(pairlist,symlist)




		

def uniqueBps(N):
	combos = list(it.combinations([i for i in range(N)],N/2))
	permlist = []
	for i in range(len(combos)/2):
		permlist.append(combos[i] + combos[len(combos)-1-i])
	return combos,permlist

def swapSpins(bp,i,j):
	i1 = bp.index(i)
	i2 = bp.index(j)
	bp[i1] = j
	bp[i2] = i
	return bp

def standardForm(bp,permlist):
	A = bp[0:len(bp)//2]
	B = bp[len(bp)//2:]
	A.sort()
	B.sort()
	if tuple(A+B) in permlist:
		return A+B
	else:
		return B+A


def equivalenceList(pairlist,symlist):
	s1 = pairlist[0][0]
	elist = [s1]
	for p in pairlist:
		if s1 in p:
			elist.append(p[(p.index(s1) + 1)%2])
	print(elist)







if __name__ == "__main__":
    main()











