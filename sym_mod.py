#!/usr/bin/env python
# coding: utf-8



import tfim
import numpy as np
from itertools import groupby
import symmetry_types

def grouped_configs(basis, Jij):
    '''function returns arrays of integers that have the same energy when converted to 
    spin configurations. excludes integers that are spin inversions of each other'''
    Energies = -1 * (tfim.JZZ_SK_ME(basis,Jij))
    sorted_indices = np.argsort(Energies)
    sorted_energies = Energies[sorted_indices]
    
    split_list = []                                #Breaks up energy list by value
    for key, group in groupby(sorted_energies):
        split_list.append(list(group))
                                            
    num_list = []
    for x in range(0, len(split_list)):           
        index = 0
        for y in range(0, x+1):
            index = index + len(split_list[y]) 
        start = (index - len(split_list[x]))
        entry = sorted_indices[start:index]
        entry.sort()
        #entry = entry[:len(entry)/2]
        num_list.append(entry)
        #num_list is a list of lists that contain indices with the same energy, excluding spin inversions
    return num_list



def list_transformations(grouped_configs, max_energy, N, basis):
    """function takes in grouped_configs, which is a list of configuration integers grouped by energy level, and
    returns a list of all transformations that occur in each energy level 
    
    when only finding transformations for a few energy levels, the ground state corresponds to max_energy=1, using
    max_energy=0 will go through all energy levels"""
    
    transformation_list = []
    if max_energy == 0:
        max_energy = len(grouped_configs)
    for entry in grouped_configs[:max_energy]:
        single_list = []
         
        if len(entry)>1: 
            index = 0
            for num in entry[:-1]:
                index = index + 1
                num = basis.state(num)
                
                for num2 in entry[index:]:
                    
                    num2 = basis.state(num2)
                    diff = num - num2
                    
                    nonzero = np.nonzero(diff)[0]
                    if nonzero.size > N/2:         #This is the piece that will throw out hamming distances greater than N/2
                        continue
                        
                    else:
                        substate1 = basis.index(num[nonzero])
                        substate2 = basis.index(num2[nonzero])

                        transformation = [nonzero, substate1, substate2]
                        single_list.append(transformation)

                        
                    #single_list.append(transformation)
        transformation_list.append(single_list)       
        #a list with nonzero array, then the integers from the binary number formed by applying the nonzero array
    return transformation_list


def symmetry(transformation_list, grouped_configs, basis):
    """function takes in a list of all transformations grouped by energy level and grouped_configs and returns a
    list of symmetries without duplicates"""
    
    symm_list = []
    level_index = 0
    for energy_level in transformation_list:
        if len(energy_level) == 0:
            continue
        level_index = level_index + 1
        transformation_index = -1
        for transformation in energy_level:
            transformation_index = transformation_index + 1
            switched_index = transformation [0]
            int1 = transformation[1]
            int2 = transformation[2] #These are integers from the binary number taken out by switched_index
            
            is_a_symmetry = True
            
            for energy_level2 in grouped_configs[:]: #used to start at level_index I think
                if  not is_a_symmetry: 
                    break
                
                for configuration in energy_level2:   #configuration is the integer
                    binary = basis.state(configuration) #convert to binary
                    switched = basis.index(binary[switched_index])  #convert flipped spins to integers like config1/2
                    if switched == int1 or switched == int2:
                        #now flip the spins in those spots and see if the new configuration is contained in
                        for value in switched_index:
                            binary[value] = (binary[value]-1)*-1
                        binary = basis.index(binary)  #now, since we have done the transformation, if this value is not in the energy, we can throw out the symmetry
                        if binary in energy_level2:
                            break
                        else:
                            is_a_symmetry = False
                            break 
                    else:
                        continue
            if is_a_symmetry:
                repeat = False
                for trans in symm_list:
                        if np.array_equal(transformation[0], trans[0]) and transformation[1]==trans[1]:
                            repeat = True
                            break
                if repeat == False:
                    symm_list.append(transformation)
                        
    return symm_list



def list_sorter(transformations, N):
    '''Takes a list of transformations or symmetries and sorts them by hamming distance
    Returns a list of symmetries sorted into lists by hamming distance'''
    sorted_list = []
    for length in range(2, N):
        one_length = []
        for transformation in transformations:
            array = transformation[0]
            if array.size == length:
                    one_length.append(transformation)
            else:
                continue
        sorted_list.append(one_length)
    return sorted_list



def remove_combo(sorted_list):
    '''Function removes transformations of Hamming distance 4 that are combinations of hamming distance 2 transformations
    Eventually, I might modify it to remove combinations that are longer though that might not become a big enough issue
    Input must be a list with Hamming distance 2 at index 0 and distance 4 at index 1'''
    
    length_2 = sorted_list[0]
    length_4 = sorted_list[1]
    transform_index = 0
    for transformation in length_2:
        transform_index += 1
        array = transformation[0]
        for transformation2 in length_2[transform_index:]:
            array2 = transformation2[0]
            if (array == array2).any() or array[0]==array2[1] or array[1]==array2[0]:
                continue
            else:
                combo = np.concatenate((array, array2))
                sorted_indices = np.argsort(combo)
                sorted_combo = combo[sorted_indices]
    
                for arr in transformation[1:]:
            
                    arr = basis.state(arr)[-2:]
                    arr2 = basis.state(transformation2[1])[-2:]
                    binary = np.concatenate((arr, arr2))
                    sorted_indices = np.argsort(combo)
                    binary = binary[sorted_indices]
                    #print 'sorted binary: ', binary
                    int1 = basis.index(binary)
                    #print 'int1: ', int1
                    
                    t_combo_index = -1
                    for t_combo in sorted_list[1]:
                        t_combo_index += 1
                        #print 't_combo: ',t_combo
                        if (t_combo[0] == sorted_combo).all() and (int1==t_combo[1] or int1==t_combo[2]):
                            del sorted_list[1][t_combo_index]
                            break
                        else:
                            continue
    return sorted_list          



def symmetry_sorter(symmetries, sortedsym):
    '''Sorts the symmetries I have found into swap, anti-swap, and other'''
    swap_symmetries = []
    anti_swap = []
    other = []
    
    swap_sym = sortedsym[0]
    anti_swap_sym = sortedsym[1]
    
    sorted_swap = []
    sorted_antiswap = []
 
    for tup in swap_sym:
        tup = sorted(tup)
        sorted_swap.append(tuple(tup))
        
    
    for tup2 in anti_swap_sym:
        tup2 = sorted(tup2)
        sorted_antiswap.append(tuple(tup2))
        
    for length in symmetries:
        for transformation in length:
            array = tuple(transformation[0])
            if array in sorted_swap:
                swap_symmetries.append(transformation)
            elif array in sorted_antiswap:
                anti_swap.append(transformation)
            else:
                other.append(transformation)
    if len(swap_symmetries) != 0:
        print 'swap_symmetries: ', swap_symmetries
    if len(anti_swap) != 0:
        print 'anti_swap: ', anti_swap
    if len(other) != 0:
        print 'other: ', other
    return swap_symmetries, anti_swap, other






