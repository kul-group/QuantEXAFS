#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from ase.io import read, write
import larch
from larch.io import read_athena
from larch_plugins.xafs import feffdat
from larch import Group, isNamedClass
from larch.xafs import feffit, TransformGroup, FeffitDataSet, feffit_report, feffit_transform, pre_edge
from larch.fitting import guess, group2params, param_group
from larch import Group, Parameter
#from larch_plugins.io<http://larch_plugins.io> import read_ascii
from larch_plugins.xafs import FeffPathGroup
from larch_plugins.xafs.feffit import (TransformGroup, FeffitDataSet,
                                       feffit, feffit_report)
from lmfit import Parameters, Parameter as param, Minimizer
from larch.xafs import feffrunner, feffpath, feff6l, feff8l
from larch.xafs import autobk
import matplotlib.pyplot as plt
from wxmplot.interactive import plot
import wx
from larch.wxlib.xafsplots import plot_chifit
import glob
import re
import sys, os
from ase import io


def read_experimental_data(filename,verbose=False, plot_expt = False):
    """
    Function to read in the experimental *.prj file
    """
    project = read_athena(filename)
    expt_data = {}
    renaming_group_map = {}
    group_count = 0
    for name, data in project._athena_groups.items():
        
        autobk(data.energy, data.mu, group=data, rbkg=0.85, kweight=3)
        #print(dir(data))

        expt_data[name] = data
        only_data= {k:v for k,v in expt_data.items() if k.startswith('PdMgO_')}

        #sys.exit()
        #autobk(data.energy, data.mu, group=data, rbkg=0.80, kweight=3)

        #print(name)

    #print(data)
    #print(data)
    #sys.exit()

    if verbose:
        for attr in dir(data):
            print(attr, type(getattr(data, attr)))

    if plot_expt:
        # Rachita - pease fix.
        plt.plot(data.k, data.chi*data.k**2, label='$\chi$')
        plt.xlabel(r'$k\, ({\rm\AA}){-1}$')
        plt.ylabel(r'$k^2\chi, ({\rm\AA})^{-2}$')
        plt.legend()

    return expt_data, only_data

def atoms2cluster(atoms,absorbing_atom='Pd',distance_cutoff=8.0):
    #atoms = io.read(filename)
    absorbing_atom_index= [a.index for a in atoms if a.symbol == absorbing_atom][0]
    atoms = atoms.repeat([2,2,1])
    atoms_cluster = atoms[[a.index for a in atoms if atoms.get_distance(absorbing_atom_index, a.index, mic=True) < distance_cutoff]]
    absorbing_atom_index = [a.index for a in atoms_cluster if a.symbol == absorbing_atom][0]
    #view(atoms_cluster)
    return atoms_cluster, absorbing_atom_index

def cluster2feffinp(atoms_super, file_name = 'feff.inp', distance_cutoff = 8.0, absorbing_atom=''):
    title='Pd_MgO'
    core = 'Pd1'
    edge = 'K'
    HOLE = edge2hole(edge)
    edge_energy = 24350

    a, b, c = atoms_super.cell.lengths()
    alpha, beta, gamma = atoms_super.cell.angles()
    absorbing_atom_index = [a.index for a in atoms_super if a.symbol == absorbing_atom][0]

    ipot = []
    distance = []
    positions = []
    atom_symbol_index_list = []
    for atom in atoms_super:
        distances = atoms_super.get_distance(absorbing_atom_index, atom.index, mic= True)
        tag = np.unique(atoms_super.get_chemical_symbols())
        Z = np.unique(atoms_super.get_atomic_numbers())
        elements_list = atoms_super.get_chemical_symbols()
        distance.append(distances)
        positions.append(atoms_super.get_distance(absorbing_atom_index, atom.index, vector=True, mic=True))
        if atom.symbol == absorbing_atom:
            index = 1
            symbol = atom.symbol
            string_rep = str(symbol) + str(index)
            atom_symbol_index_list.append(string_rep)
        else:
            index = atom.index +2
            symbol = atom.symbol
            string_rep = str(symbol) + str(index)
            atom_symbol_index_list.append(string_rep)
        final_str = atom_symbol_index_list


    for i in elements_list:
        if i == 'O':
            i_pot = 2
        else:
            if i == absorbing_atom:
                i_pot = 0
            else:
                i_pot = 1
        ipot.append(i_pot)
    
    file_text = f"""
 * title = name:     {title}
 * title = formula:  {title}
 * space = p 42/m n m
 * a     =   {a}    b    =   {b}    c     =   {c}
 * alpha =  {alpha}    beta =  {beta}    gamma =  {gamma}
 * rmax  =   {distance_cutoff}    core  = {core}
 * polarization = 0  0  0
 * shift =   0  0  0
 * atoms
 * # el.        x           y           z       tag
 *   {absorbing_atom}     0.00000     0.00000     0.00000     {core}
 *  
 * --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--
 * --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--

 TITLE name:     {title}
 
 HOLE      {HOLE}   1.0  *  FYI: ({core} {edge} edge @ {edge_energy} eV, 2nd number is S0^2)
 *         mphase,mpath,mfeff,mchi
 CONTROL   1      1     1     1  
 PRINT     1      0     0     0  
 RMAX      5.0
  * POLARIZATION  0   0   0
 POTENTIALS
  * ipot   Z      tag
     0     {Z[2]}     {tag[2]}
     1     {Z[1]}     {tag[0]}
     2     {Z[0]}     {tag[1]}
 ATOMS                  * this list contains {len(atoms_super)} atoms
 *  x              y              z        ipot      tag      distance
 """

    for atom, t, i, d in zip(atoms_super, final_str, ipot, distance):
        file_text +=f"   {positions[atom.index][0]:.5f}       {positions[atom.index][1]:.5f}       {positions[atom.index][2]:.5f}     {i}       {t}       {d:.5f}\n"
    file_text += "END\n"
    with open (file_name, 'w') as file:
        file.write(file_text)

    print('Finished atoms2feff_inp')
def edge2hole(edge):
    '''
    A quick way to convert the string name of an edge to a hole value for feff6
    use in conjunction with cluster2feffinp

    Parameters
    ----------
    edge : str
        the edge that feff will calcualte scattering parameters at

    Returns
    -------
    hole : int
        interger value used by feff 6 to represent which excitation event.

    '''
    edges = ['K', 
             'L1', 'L2', 'L3', 
             'M1', 'M2', 'M3', 'M4', 'M5',
             'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7',
             'O1', 'O2', 'O3', 'O4', 'O5',
             'P1', 'P2', 'P3']
    
    hole = edges.index(edge) + 1
    return hole
    HOLE = edge2hole(edge = 'K')
def read_paths():
    paths = []
    pars = param_group(
              del_e0 = param(0.7, vary=True),
              sig2_1 = param(0.002, vary=True),
              sig2_2 = param(0.003, vary=True),
              sig2_3 = param(0.004, vary=True),
              sig2_4 = param(0.004, vary=True),
              sig2_5 = param(0.005, vary=True),
              sig2_6 = param(0.004, vary=True),
              del_r = param(0.0, vary=True))
    file = open('files.dat')
    file_list = file.readlines()[9:]
    vals = []
    i=1
    for f in file_list:
        vals.append(f.split())
        paths[i]= feffdat.feffpath(f, sigma2='sig2', degen = 'deg', deltar='del_r*reff')


        paths[i]= feffdat.feffpath(f, sigma2='sig2_1', degen = 'deg', deltar='del_r*reff')
        paths[i]= feffdat.feffpath(f, sigma2='sig2_2', degen = 'deg', deltar='del_r*reff')


        i+=1
    print (paths)
    return (paths)

#reading values for path parameters from files.dat and passing them to paths in read_paths
"""filename = []
filename = glob.glob('*.dat', recursive = True)
for f in filename:
    def get_pars(filename):
        with open('files.dat') as file:
            line = [l for l in file if filename in l][0]
            sig2, amp_ratio, deg, nlegs, r_eff = list(map(float, line.split()[1:]))
            pathargs = sig2, amp_ratio, deg, nlegs, r_eff
        return pathargs"""

#path1 = feffdat.feffpath ('feff0001.dat', degen = 'n3*6', e0 = 'del_e0', sigma2 = 'sig2_1', deltar = 'del_r*reff')
#formatting the paths for reading parameters and making paths readable for the transform function
def trans_func(paths):
    trans = TransformGroup(kmin =3.0, kmax=11.5, kw=3, dk=5, window = 'hanning', rmin=1.2, rmax=5)
    for path in paths:
        dset = FeffitDataSet(data = read_data.data, pathlist=[paths], transform=trans)
        out = feffit(read_paths.pars, dset)
    return (out)


def report(out):
    report = feffit_report(out)
    file_name = 'report.txt'
    with open (file_name, 'w') as file:
        file.write(report)
    return

def plot_fits(dset):
    app = wx.App()
    plot_chifit(dset, rmax = 5.0, show_mag=True, show_real=True)
    return




#using atoms code to create feff.inp file and running feff calculations
def atoms2feff():
    feff_inp = []
    feff_list = glob.glob('*.cif', recursive = False)
    for f in feff_list:
        feff_inp[f] = atoms2feff.write_input()
    return ()

"""reading the Athena prj file and printing the contents of the project file and also the contents
each group in the file and doing background substraction using autobk module"""

