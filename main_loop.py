#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 00:33:55 2021

@author: rachita
"""
from larch.io import read_athena, write_group, create_athena
from larch.xafs import autobk
import os, sys
import numpy as np
from try_1 import exafs_fitting, atoms2report
import matplotlib.pyplot as plt



def read_experimental_data(filename,verbose=False, plot_expt = False):
        """
        Function to read in the experimental *.prj file
        """
        project = read_athena(filename)
        print ('beginning to read the data')
        expt_data = {}
    
        renaming_group_map = {}
        group_count = 0
    
        for name, data in project._athena_groups.items():
            
            autobk(data.energy, data.mu, group=data, rbkg=0.85, kweight=3)
            expt_data[name] = data
            only_data= {k:v for k,v in expt_data.items() if k.startswith('PdMgO_')}
            if plot_expt:
            # Rachita - pease fix.
                plt.plot(data.k, data.chi*data.k**2, label='$\chi$')
                plt.xlabel(r'$k\, ({\rm\AA}){-1}$')
                plt.ylabel(r'$k^2\chi, ({\rm\AA})^{-2}$')
                plt.legend()
    #################################################
        #only_data= {k:v for k,v in expt_data.items() if k.startswith('Ga_')}           
    
    
        if verbose:
            for attr in dir(data):
                print(attr, type(getattr(data, attr)))
            
        
        return expt_data, only_data
    
expt_data, only_data = read_experimental_data('/Users/rachita/Box/Larch_part2/Pd_atom_dis/Insitu_EXAFS_temp/PdMgO.prj')
print(only_data)

wd2 = os.getcwd()
i = 0
for k,v in only_data.items():#range (len(only_data)):
    path = '{}'.format(k)
    print (v)
    if not os.path.exists(path):
        os.mkdir(path)
        os.chdir(path)
        exafs_fitting(v)
        
    os.chdir(wd2)
########################################################
      