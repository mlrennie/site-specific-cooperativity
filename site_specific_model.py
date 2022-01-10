__description__ = \
"""
site_specific_model.py

Class for site-specific model description based on a hierarchical parameterization.

Warning: computation speed is extremely dependent on number of sites. For models
with more than 8 sites computation maybe longer than one minute to generate the 
binding polynomial. For optimum speed apply symmetry and simplifications first 
then generate equations.

Note that once simplifications are applied they cannot be undone within the same model.
"""
__author__ = "Martin L. Rennie"
__date__ = "2021-02-24"

from site_specific_utilities import *

class SiteSpecificModel:
    '''
    Class for site-specific models description based on a hierarchical parameterization.
    '''
    
    def __init__(self, n_sites, homotropic=False, sym='None', n_asu=None):
        '''
        n_sites : int, number of binding sites
        sym : string, symmetry of the system (currently only cyclic)
        homotropic : bool, set to true if all sites bind the same ligand
        '''
        if(type(n_sites) is not int or n_sites<1):
            raise ValueError("Number of sites number be a positive integer.")
            
        self.n_sites = n_sites
        self.states = generate_states(n_sites)
        self.grouped_states = group_states(self.states)
        self.equiv_states = []

        self.params = {}
        
        self.homotropic = homotropic
        self.L = [None]*self.n_sites
        self.make_ligand()
        
        if(n_asu == None):
            self.n_asu = self.n_sites
        else:
            self.n_asu = n_asu
        self.sym = sym.lower()
        self.apply_sym()
        
        self.make_parameters()
        
    def make_ligand(self):
        '''
        Generate a list of sympy variables for each ligand.
        '''
        for i in range(0,self.n_sites):
            if(self.homotropic):
                self.L[i] = symbols('[L]')      
            else:
                self.L[i] = symbols('[L_{{{}}}]'.format(i+1))
    
    def make_parameters(self):
        '''
        Generate a dictionary of model parameters.
        '''
        # try improve speed by converting to binary
        for s in self.states[1:]:
        
            # determine the number of occupied sites for the state
            # which will take on values of 0 or 1 in the product
            n_occ = s.count('1')
        
            # find the indices of the unoccupied sites
            indices = [m.start() for m in re.finditer('0',s)]
        
            # find all the sub-states for the given number of occupied sites
            s_sub = [''.join(seq) for seq in 
                      itertools.product('01', repeat=n_occ)]
            s_sub.reverse()
        
            # iterate through the sub-states and use these to populate the parameters
            term = 1
            for sub in s_sub:
                # add in the unoccupied sites
                for i,index in enumerate(indices):
                    sub = sub[:index]+'0'+sub[index:]
                sub_n_occ = sub.count('1')
                # merge any equivalenced states
                if(self.equiv_states!=[]):
                    for j,equiv in enumerate(self.equiv_states[sub_n_occ]):
                        if(sub in equiv):
                            sub = self.nr_states[sub_n_occ][j]
                # make sympy variables for the parameters
                if(sub not in self.params):
                    if(sub.count('1')==1):
                        #self.params.update({sub:symbols('K_{{{}}}'.format(int(sub,2)))})
                        self.params.update({sub:symbols('K_{{{}}}'.format(sub))})
                    elif(sub.count('1')>1):
                        #self.params.update({sub:symbols(r'\alpha_{{{}}}'.format(int(sub,2)))})
                        self.params.update({sub:symbols(r'\alpha_{{{}}}'.format(sub))})
                        
    def make_params_list(self):
        '''
        Make a list of the equilibrium constants of the model.
        '''
        # initialize params_list with non-parameter values
        self.params_list = [0,1]
        for param in self.params.values():
            if(param not in self.params_list):
                self.params_list.append(param)
        # remove non-parameter values from list
        self.params_list = self.params_list[2:]
          
### Binding polynomial and associated equations ###       
    
    def generate_bp_terms(self):
        '''
        Generate a dictionary of sympy expressions for each state.
        '''
        ##slow, speed up by allowing approximations first ##
        check_states(self.states)
        
        if(self.L==[None]*self.n_sites):
            self.make_ligand()
        
        try:
            self.terms
        except AttributeError:
            self.terms = {'0'*self.n_sites:1}
        
        for s in self.states[1:]:
        
            # determine the number of occupied sites for the state
            # which will take on values of 0 or 1 in the product
            n_occ = s.count('1')
        
            # find the indices of the unoccupied sites
            indices = [m.start() for m in re.finditer('0',s)]
        
            # find all the sub-states for the given number of occupied sites
            s_sub = [''.join(seq) for seq in 
                      itertools.product('01', repeat=n_occ)]
            s_sub.reverse()
        
            # iterate through the sub-states and put these in the product
            term = 1
            for sub in s_sub:
                # add in the unoccupied sites
                for i,index in enumerate(indices):
                    sub = sub[:index]+'0'+sub[index:]
                sub_n_occ = sub.count('1')
                # merge any equivalenced states
                if(self.equiv_states!=[]):
                    for j,equiv in enumerate(self.equiv_states[sub_n_occ]):
                        if(sub in equiv):
                            sub = self.nr_states[sub_n_occ][j]
                # compute the term based on f(b) 
                if(sub.count('1')==1):
                    i = sub.find('1')
                    term = Mul(Mul(term,self.params[sub]),self.L[i]) # without simplifications this step is slow
                elif(sub.count('1')>1):
                    term = Mul(term,self.params[sub]) # without simplifications this step is slow
            
            self.terms.update({s:term})
    
    def generate_bp(self):
        '''
        Generate a sympy expression for the binding polynomial.
        '''
        
        try:
            self.terms
        except AttributeError:
            self.generate_bp_terms()
        
        self.bp = sum(self.terms.values())
    
    def generate_mass_balance(self):
        '''
        Generate sympy expressions for the mass equations.
        Based on Eqns 26 from Vega et al, 2015, Methods.
        '''
        
        # need binding polynomial to compute mass balance equations
        try:
            self.bp
        except AttributeError:
            self.generate_bp()
                
        self.mb_eqns = {}
        M0 = symbols('[M_{{{}}}]'.format(self.states[0]))
        
        if(self.homotropic):
            lig = self.L[0]
            ligtot = 'L_t'
            self.mb_eqns.update({'L_t':diff(self.bp,lig)*lig*M0})
        else:
            for i,lig in enumerate(self.L):
                ligtot = 'L_{}t'.format(i+1)
                self.mb_eqns.update({ligtot:diff(self.bp,lig)*lig*M0})
            
        self.mb_eqns.update({'M_t':self.bp*M0})
    
    def generate_fractional_saturation(self):
        '''
        Generate sympy expressions for the fractional saturation 
        of each ligand
        '''
        
        # need binding polynomial to compute mass balance equations
        try:
            self.bp
        except AttributeError:
            self.generate_bp()
            
        ####### TO COMPLETE #####
    
    def reset_eqns(self):
        '''
        Recompute the binding polynomial, mass balance, and heat equations.
        '''
        
        # if terms, binding polynomial, mass balance, or heat equation have been generated 
        # re-compute them
        try:
            self.terms
            self.generate_bp_terms()
        except AttributeError:
            pass
        
        try:
            self.bp
            self.generate_bp()
        except AttributeError:
            pass
            
        try:
            self.mb_eqns
            self.generate_mass_balance()
        except AttributeError:
            pass
        
### Simplifications ###
    
    def apply_sym(self):
        '''
        Generate a list of equivalenced states based on a specified symmetry.
        Only single site symmetry currently implemented.
        
        sym : string, 'cyclic' for cyclic symmetry
        n_asu : positive integer, number of asymmetric units making up the host.
        '''
        n_asu_sites = self.n_sites/self.n_asu
        if(not n_asu_sites.is_integer() or n_asu_sites<0):
            raise ValueError("Number of sites must be divisible by "
                             "the number of asymmetric units") 
        
        if(self.sym=='cyclic'):
            self.equiv_states,self.nr_states = apply_cyclic_sym(self.grouped_states,
                                                                self.n_sites,self.n_asu)
            # force ligands for equivalent sites to be identical
            for i in range(int(n_asu_sites)):
                for j in range(1,self.n_asu):
                    self.L[i*self.n_asu+j] = self.L[i*self.n_asu]
        elif(self.sym=='dihedral'):
            if(self.n_sites%2!=0):
                raise ValueError("Even number of sites required for dihedral symmetry.")
            if(self.n_asu==self.n_sites):
                self.equiv_states,self.nr_states = apply_dihedral_sym(self.grouped_states,
                                                                      self.n_asu)
                self.homotropic = True
                self.make_ligand()
            else:
                raise ValueError("Multiple sites per asymmetric unit not currently implemented"
                                 " for dihedral symmetry.")
        elif(self.sym=='none'):
            self.equiv_states = []
        else:
            raise ValueError("Symmetry not recognised or has not been implemented yet.")
        
        # if terms, binding polynomial, mass balance, or heat equation have been generated 
        # re-compute them
        self.reset_eqns()

    def nearest_neighbors(self):
        '''
        Apply strict nearest-neighbor approximations.
        
        This assumes adjacent binding sites in the structure are adjacent
        in the binary notation.
        '''
        for param in self.params:
            if(param.count('1')>2):
                self.params[param] = 1
            elif(param.count('1')==2 and param.count('11')!=1):
                self.params[param] = 1
    
        # if terms, binding polynomial, or mass balance equations have been generated re-compute them
        self.reset_eqns()
            
        # if enthalpies have been produced constrain them and reproduce the heat equation
        ##### TO DO ####
    
    def na_nearest_neighbors(self):
        '''
        Apply non-additive nearest-neighbor approximations.
        This assumes adjacent binding sites in the structure are adjacent
        in the binary string.
        '''
        for param in self.params:
            if(param.count('1')>3):
                self.params[param] = 1
            elif(param.count('1')==3 and param.count('111')!=1):
                self.params[param] = 1
            elif(param.count('1')==2 and param.count('11')!=1):
                self.params[param] = 1
    
        # if terms, binding polynomial, or mass balance equations have been generated re-compute them
        self.reset_eqns()
        
    def set_equivalent_affinities(self,sites):
        '''
        Make a group of sites equivalent in terms of binding affinity.
        
        sites : list of positive integers, these are the binding site
                positions to be equivalenced as ordered in the binary notation.
                Note that the first position is 1, not 0 as per Python lists.
        '''
        # order from smallest to largest
        sites.sort()
        # set the reference 
        s_ref = '0'*(sites[0]-1)+'1'+'0'*(self.n_sites-sites[0])
        # equivalence association constants
        for site in sites[1:]:
            s = '0'*(site-1)+'1'+'0'*(self.n_sites-site)
            self.params[s] = self.params[s_ref]
            
    def set_independent(self,s_list):
        '''
        Simplify the model by assuming no linkage between specified sets of sites.
        
        s_list : list of binary strings representing the linkage constants to be fixed
        '''
        for param in self.params:
            if(param in s_list):
                self.params[param] = 1
                
        self.reset_eqns()
            
    def set_block(self,s_list):
        '''
        Simplify the model by assuming extreme negative linkage between specified sets of sites.
        
        s_list : list of binary strings representing the linkage constants to be fixed
        '''
        for param in self.params:
            if(param in s_list):
                self.params[param] = 0
        
        self.reset_eqns()
        