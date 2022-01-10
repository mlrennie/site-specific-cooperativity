__description__ = \
"""
site_specific_utilities.py

associated functions to assist site_specific_model.py

"""
__author__ = "Martin L. Rennie"
__date__ = "2021-02-24"

import itertools
import re
from matplotlib import pylab as plt
import numpy as np
from symengine import *; from sympy import simplify


def generate_states(n_sites):
    '''
    Generate a list of binary strings of a given string length.
    '''
    to_return = [''.join(seq) for seq in 
                  itertools.product('01', repeat=n_sites)]
    return to_return
    
def check_states(states):
    '''
    Ensure list contains only binary vectors.
    '''
    if(not type(states) is list):
        raise Exception("Input should be a list")
    
    n_sites = len(states[0])
    reg = re.compile(r'[^0-1]')
    for state in states:
        if(bool(reg.search(state))):
            raise Exception("Error reading state: ",
                            state,
                            "\nAll states must be binary vectors.")
        if(len(state)!=n_sites):
            raise Exception("Error reading state: ",
                            state,
                            "\nAll binary vectors must be of equal length.")

def group_states(states):
    '''
    Group a list of binary strings by the number of 1's
    '''
    # ensure correct list semantics
    check_states(states)
    
    n_sites = len(states[0])
    grouped_states = [[] for i in range(n_sites+1)] # initialize list
    
    for state in states:                            # count number
        n_ligands = state.count('1')                # of occupied sites
        grouped_states[n_ligands].append(state)
    
    [i.reverse() for i in grouped_states]
    
    return grouped_states

def cycle_elements(l, n1, n2):
    '''
    Cycle the elements of a list, l, from right to left by n1 spaces, in groups of n2
    '''
    out=''
    for i in range(int(len(l)/n2)):
        out = out+l[n1+i*n2:(i+1)*n2]+l[i*n2:n1+i*n2]
    return out

def cycle_elements_dihedral(l,n):
	'''Cycle the elements of a list, l, from right to left by n spaces. Comparible with dihedral format'''
	
	half = (len(l))//2
		
	return l[:half][n:] + l[:half][:n] + l[half:][n:] + l[half:][:n]

def step_down(state):
    '''
    Decompose a binary string into a list of strings with one less '1'
    '''
    tmp = []
    for i in range(len(state)):
        if(state[i]=='1'):
            tmp.append(state[:i]+'0'+state[i+1:])
    return tmp

### simplify sympy binding polynomial using symmetry/approximations
def apply_cyclic_sym(grouped_states,n_sites,n_asu):
    '''
    Apply cyclic symmetry to binary vectors, returning a list 
    of equivalent states.
    '''
    
    equiv_states = [[] for i in range(n_sites+1)]
    nr_states = [[] for i in range(n_sites+1)]
    
    for i,states in enumerate(grouped_states):
        
        # create a copy as the base for a new list
        dup_states = states.copy()
        
        # cycle through states with the same number of ligands
        # bound until all states have been accounted for
        while(dup_states != []):
            k = 1
            perm_state = []
            
            equiv_states[i].append([dup_states[0]])
            
            # cycle through symmetry-equivalent states until returning to the start
            while(dup_states[0] != perm_state):
                
                perm_state = cycle_elements(dup_states[0],k,n_asu)
                
                # if the permutation is within the list of states account
                # for it by removing it and placing it in the equivalenced states
                if(perm_state in dup_states[1:]):
                    # add the permutation state to the list of equivanced states
                    equiv_states[i][-1].append(perm_state)
                    dup_states.remove(perm_state)
                k = k + 1
            
            nr_states[i].append(dup_states[0])
            equiv_states[i][-1] = set(equiv_states[i][-1]) # use sets for speed
            del(dup_states[0])
    
    return equiv_states,nr_states

def apply_dihedral_sym(grouped_states,n_asu):
    '''
    Apply cyclic symmetry to binary vectors, returning a list 
    of equivalent states.
    '''
    
    equiv_states = [[] for i in range(n_asu+1)]
    nr_states = [[] for i in range(n_asu+1)]
    
    for i,states in enumerate(grouped_states):
        
        # create a copy as the base for a new list
        dup_states = states.copy()
        
        # cycle through states with the same number of ligands
        # bound until all states have been accounted for
        while(dup_states != []):
            #k = 1
            perm_state = []
            
            equiv_states[i].append([dup_states[0]])
            
            # cycle through symmetry-equivalent states until returning to the start
            n_asu_half = n_asu//2
            for k in range(1,n_asu+1):
                
                if(k < n_asu_half):
                    perm_state = cycle_elements_dihedral(dup_states[0],k)
                else:
                    perm_state = cycle_elements_dihedral(dup_states[0][::-1],k-n_asu_half)
                # if the permutation is within the list of states account
                # for it by removing it and placing it in the equivalenced states
                if(perm_state in dup_states[1:]):
                    # add the permutation state to the list of equivanced states
                    equiv_states[i][-1].append(perm_state)
                    dup_states.remove(perm_state)
              
            nr_states[i].append(dup_states[0])
            equiv_states[i][-1] = set(equiv_states[i][-1]) # use sets for speed
            del(dup_states[0])
        
    return equiv_states,nr_states

def plot_schematic(grouped_states,n_sites):
    '''
    Plot the stepwise equilibrium scheme.
    Returns fig and ax from subplots.
    '''
    # plot the grouped states
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))   
        
    max_sub = 0
    group_tmp1, group_tmp2 = ([],[])
    offset, offset2 = (0.25,0.03)   # emperical offsets
    for x,group in enumerate(grouped_states):
        n_sub = len(group)
        if(n_sub>max_sub):
            max_sub = len(group)
            
        for j,state in enumerate(group):
            y = n_sub/2-j-0.5
            ax.text(x,y,state,fontsize=10,
                    horizontalalignment='center',
                    verticalalignment='center',
                    bbox=dict(facecolor='none', edgecolor='black', 
                              boxstyle='round4', fc='white'))
            
            # plot equilibrium arrows
            group_tmp1.append([state,x,y])
            down_states = step_down(state)
            
            for group_tmp in group_tmp2:
                if(group_tmp[0] in down_states):
                    
                    # empirical offset to 'scale' arrows
                    if(group_tmp[2]<y):
                        offset3 = 0.1
                    elif(group_tmp[2]==y):
                        offset3 = 0.
                    else:
                        offset3 = -0.1
                    
                    # forward arrow
                    ax.arrow(group_tmp[1]+offset,group_tmp[2]+offset2+offset3,
                             x-group_tmp[1]-offset*2,
                             y-group_tmp[2]-2*offset3,
                             head_width=0.05,head_length=0.05,
                             length_includes_head=True,shape='right',
                             fc='k',linewidth=0.2)
                    # reverse arrow
                    ax.arrow(x-offset,y-offset2-offset3,
                             group_tmp[1]-x+offset*2,
                             group_tmp[2]-y+2*offset3,
                             head_width=0.05,head_length=0.05,
                             length_includes_head=True,shape='right',
                             fc='k',linewidth=0.2)
        group_tmp2 = group_tmp1
    
    ax.set_xlim(-0.5,n_sites+0.5)
    ax.set_ylim(-max_sub/2.,max_sub/2.)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Protein States')
    
    fig.set_size_inches(3+n_sites, max_sub)
    
    return (fig, ax)

def sympy_to_function(self, sympy_vars, func):
		'''
		Convert sympy binding polynomial expression to a function.
		
		sympy_vars : ordered list of sympy variables for the function
		func	   : sympy expression to be converted to a function
		'''
		return lambdify(sympy_vars,func) 