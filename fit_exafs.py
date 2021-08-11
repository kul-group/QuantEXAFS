#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
import sys
import os
from ase.visualize import view
from utils import read_experimental_data, atoms2cluster, cluster2feffinp
from ase.db import connect
from numpy import gradient, ndarray, diff, where, arange, argmin
import matplotlib.pyplot as plt
from larch.wxlib import plotlabels as plab

def exafs_fitting (data1):
    import os
    from ase.io import read, write
    from ase.db import connect
    import traceback
    db = connect('/Users/rachita/Box/Larch_part2/Pd_atom_dis/Insitu_EXAFS_temp/pd_mgo.db')
    wd = os.getcwd()
    print ('Reading the ase.db')
    #from exafs_fit import atoms2report
    c= data1 
    
    for row in db.select(metal='Pd'):
        facet = row.facet
        sub = row.sub
        ads = row.ads
        config = row.config 
        name = f"{facet}_{sub}_{ads}_{config}"
        print(name)
        if not os.path.exists(name):
            os.mkdir(name)
            os.chdir(name)
            atoms = row.toatoms()
            if name == 's100_sub1_O2_0': 
                pass
            else:
                atoms.write(name+'.cif')
               
                try:
                    atoms2report(atoms, c)
                    #print(os.getcwd())
                except:
                    print ('Failed to read object or write to db')
            os.chdir(wd)
            
        else:
            pass   
    return 


def atoms2report(atoms, c):
    # Read Athena *.prj file
    ####### EDIT THIS TO INCLUDE FULL PATH ########
    #expt_data, only_data = read_experimental_data('/Users/rachita/Box/Larch_part2/Pd_atom_dis/Insitu_EXAFS_temp/PdMgO.prj')
    #print(only_data)
    for attr in dir(c):
        print(attr, type(getattr(c, attr)))
    
    # Read our reference atoms object and convert to feff.inp
    atoms_cluster, absorbing_atom_index = atoms2cluster(atoms,absorbing_atom='Pd',distance_cutoff=8.0)
    cluster2feffinp(atoms_cluster,absorbing_atom='Pd')
    
    #####
    #def xx
    #from lmfit import Parameter as param
    from larch import Group
    from larch.fitting import guess, group2params, param_group, param
    from larch.xafs import feffrunner, feffpath, feff6l, feff8l
    from larch.xafs import feffit, TransformGroup, FeffitDataSet, feffit_report, feffit_transform, pre_edge
    from larch_plugins.xafs import feffdat
    
    # This needs to be moved - RR
    feff6l(feffinp='./feff.inp')
    #sys.exit()
    
    # Analyze files.dat to automatically detect paths
    f = open('files.dat', "r")
    list_nlegs = []
    
    dict_sig2s = {2:'sig2_2',3:'sig2_3',4:'sig2_4',5:'sig2_5',6:'sig2_6'}
    list_pathargs = []
    for line in f.readlines():
        i_line = line.strip()
        if 'feff' in i_line:
            i_fname, i_sig2, i_amp, i_deg, i_nlegs, i_r_eff = i_line.split()
    
            var_fname = i_fname
            #if float(i_r_eff) <3.7:
            if float(i_r_eff) < 2.1: 
                var_sig2 = 'sig2_2' 
            elif float(i_r_eff) > 2.1 and float(i_r_eff)<3.2: 
                var_sig2 = 'sig2_3' 
            elif float(i_r_eff) > 3.2 and float(i_r_eff)<3.8: #and int(i_nlegs) == 2: 
                var_sig2 = 'sig2_4' 
            elif float(i_r_eff) > 3.8 and float(i_r_eff)<4.5: #and int(i_nlegs)>2: 
                var_sig2 = 'sig2_5' 
            else:
                var_sig2 = 'sig2_6'
            #var_sig2 = dict_sig2s[int(i_nlegs)]
            
            var_degen = 'deg' 
            if float(i_r_eff)<3.1:
                var_deltar = 'del_r1*reff' #del_r*%s' % i_r_eff - check
            elif float(i_r_eff)>3.1 and float(i_r_eff)<3.7:
                var_deltar = 'del_r2*reff' 
            else:
                var_deltar = 'del_r3*reff'
            #print(var_fname, var_sig2, var_degen, var_deltar)
            i_pathargs = [var_fname, var_sig2, var_degen, var_deltar]
            list_pathargs.append(i_pathargs)
    
            #print(i_fname, i_nlegs, dict_sig2s[int(i_nlegs)])
            list_nlegs.append(i_nlegs)
    f.close()
    num_paths = len(list_nlegs) #Total paths we care about for now
    print('Total paths read = %s' % num_paths)
    #return list_pathargs
    
    
   
    
    paths = []
    pars = param_group(del_e0 = param(0.7, vary=True),
              sig2_2 = param(0.002, min=0.0, max=0.1, vary=True),
              sig2_3 = param(0.003, min=0.0, max=0.1, vary=True),
              sig2_4 = param(0.004, min=0.0, max=0.1, vary=True),
              sig2_5 = param(0.005, min=0.0, max=0.1, vary=True),
              sig2_6 = param(0.005, min=0.0, max=0.1, vary=True),
              del_r1 = param(0.0, vary=True),
              del_r2 = param(0.0, vary=True),
              del_r3 = param(0.0, vary=True))
    
    for i, i_pathargs in enumerate(list_pathargs):
        var_fname, var_sig2, var_degen, var_deltar = i_pathargs
        pathargs = dict(e0='del_e0', sigma2=var_sig2, deltar=var_deltar)
        paths.append(feffdat.feffpath(var_fname,**pathargs))
        #paths[i]= feffdat.feffpath(f, sigma2=var_sig2, degen = var_degen, deltar=var_deltar)
    
    ## Fit stuff trasnform
    ###########################
    #def trans_func(paths):
    trans = TransformGroup(kmin =2.2, kmax=12.5, kweight=3, dk=1, window='hanning', rmin=1.0, rmax=5.0, fitspace='r')
    #for path in paths:
    dset = FeffitDataSet(data=c, pathlist=paths, transform=trans)
    
    out = feffit(pars, dset)
    #return (out)cd
    print(out)
    
    # Fit stuff report
    #def report(out):
    report = feffit_report(out)
    file_name = 'report.txt'
    with open (file_name, 'w') as file:
        file.write(report)
    #    return
    import numpy as np
    from scipy.interpolate import CubicSpline
    import similaritymeasures
    dataset= dset
    expt_x = dataset.data.r
    expt_y = dataset.data.chir_mag
    sim_x = dataset.model.r
    sim_y = dataset.model.chir_mag
    def compare_curves(sim_x = dataset.data.r, sim_y=dataset.data.chir_mag, expt_x=dataset.model.r, expt_y=dataset.model.chir_mag, x_range=None):
        """Returns Fréchet distance, area between curves, and RMSE
           for simulation and experimental spectra
    
           Input:
               sim_x = list of x values for simulation curve
    
               sim_y = list of y values for simulation curve
    
               expt_x = list of x values for experimental curve
    
               expt_y = list of y values for experimental curve
    
               x_range (optional) = [min_x, min_y], range of x values to use
        """
        if x_range is None: # uses full range if none supplied
            x_min = min(sim_x)
            x_max = max(sim_x)
        else:
            x_min, x_max = x_range
    
        if list(sim_x) != list(expt_x):
            cs = CubicSpline(expt_x, expt_y)
            expt_x = sim_x
            expt_y = cs(expt_x)
    
        sim = [[x, y] for x, y in zip(sim_x, sim_y) if x <= x_max and x >= x_min]
        expt = [[x, y] for x, y in zip(expt_x, expt_y) if x <= x_max and x >= x_min]
        sim, expt = np.array(sim), np.array(expt)
    
        frechet = similaritymeasures.frechet_dist(sim, expt)
        area = similaritymeasures.area_between_two_curves(sim, expt)

        mse = sum((sim[:, 1] - expt[:, 1])**2) / len(sim[:, 1])
        rmse = np.sqrt(mse)
        print ('*************************************************')
        
        return frechet, area, rmse
   
    
    compare = compare_curves(sim_x, sim_y, expt_x, expt_y, x_range=None)
    comp = f'frechet = {compare[0]}\narea = {compare[1]}\nRMSE = {compare[2]}\n'
    
    file_name = 'compare_curve.txt'
    with open (file_name, 'w') as file:
        file.write(comp)
    import seaborn as sns
    colors = sns.color_palette('muted')
    colors2 = sns.color_palette('bright')
    
            
            
    def plot_chifit(dataset, kmin=0, kmax=None, kweight=None, rmax=None,
                show_mag=True, show_real=False, show_imag=False,
                title='name', new=True, delay_draw=False, offset=0.25, win=1,
                _larch=None):
        if kweight is None:
            kweight = dataset.transform.kweight
        
        if isinstance(kweight, (list, tuple, ndarray)): kweight=kweight[0]
        data_chik  = dataset.data.chi * dataset.data.k**kweight
        model_chik = dataset.model.chi * dataset.model.k**kweight
    
        # k-weighted chi(k) in first plot window
        plt.figure()   
        plt.xlim([0,14])
        plt.plot(dataset.data.k, data_chik+offset, color='black', label='data')        
        plt.plot(dataset.model.k, model_chik+offset, color=colors[0], label='fit')
        plt.xlabel(plab.k)
        plt.ylabel(plab.chikw.format(kweight))
        plt.legend()
        plt.savefig('k-space.svg')
        # plotting the real part and the magnitude of the fit
        
        if show_mag:
            #plt.figure()
            plt.plot(dataset.data.r,  dataset.data.chir_mag+offset,
                 color='black',label='|data|')
            plt.plot(dataset.model.r, dataset.model.chir_mag+offset,
                  color=colors2[0], label='|fit|')
            plt.xlabel(plab.r)
            plt.ylabel(plab.chir.format(4))
            plt.legend()
            plt.xlim([0,5])
            plt.savefig('R-space_mag.svg')
        
        if show_real:
            plt.figure()
            plt.plot(dataset.data.r, dataset.data.chir_re+offset, color='black', label='Re|data|')
            plt.plot(dataset.model.r, dataset.model.chir_re+offset, color=colors[2], label='Re|fit|')
            plt.legend()
            plt.xlabel(plab.r)
            plt.ylabel(plab.chirre.format(3))
            plt.xlim([0,5])
            plt.savefig('R-space_re.svg')
#######################################################################################################            
            file_name = 'plot_data.txt'
            with open (file_name, 'w') as file:
                np.savetxt('plot_data', [dataset.data.r, dataset.data.chir_mag], header='expt data in r-space (real)', 
                       footer='second array is chir_mag')
            file_name = 'plot_model.txt'
            with open (file_name, 'w') as file:
                np.savetxt('plot_model', [dataset.model.r, dataset.model.chir_mag], header='model in r-space (real)', 
                       footer='second array is chir_mag')
#printing the dset values for the plots        
        if show_imag:
            plt.plot(dataset.data.r, dataset.data.chir_im+offset, linestyle='dotted', color='black', label='Im|data|')
            plt.plot(dataset.model.r, dataset.model.chir_im+offset, linestyle='dashed', color='red', label='Im|fit|')
        
        plt.figure()
        plt.plot(dataset.data.r,  dataset.data.chir_mag+offset,
                 color='black',label='|data|', alpha=0.7)
        plt.plot(dataset.model.r, dataset.model.chir_mag+offset,
                  color=colors2[0], label='|fit|', linewidth=2.0, alpha=0.8)
        plt.xlabel(plab.r)
        plt.ylabel(plab.chir.format(4))
        plt.legend()
        plt.xlim([0,5])
        plt.plot(dataset.data.r, dataset.data.chir_im+offset, color='black', label='Im|data|', alpha=0.5)
        plt.plot(dataset.model.r, dataset.model.chir_im+offset,  color=colors[2], label='Im|fit|')
        plt.legend()
        plt.savefig('R-space.svg')
    plot_chifit(dset, rmax = 5.0, show_mag=True, show_real=True)   
    
    