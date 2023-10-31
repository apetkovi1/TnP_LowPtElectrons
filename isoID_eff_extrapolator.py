import math
from turtle import setundobuffer
from matplotlib.ticker import NullFormatter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from matplotlib import rc

import numpy as np
import pandas as pn
import uproot as up
import awkward as ak
import json
import mplhep as hep

print('Libraries read in')
plt.style.use(hep.style.CMS)

path = '/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/DPS_plots/' # The output directory for the plots

############################################################

##########     Extrapolation plots

#df_eff = pn.read_csv('RelIso_bin_eff_val_HZZ.txt',sep='\t') # output from TnP tool
df_eff = pn.read_csv('/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_10/src/TnP_fits/TnP_LowPtElectrons/plots/iso_HZZ/UL_altRez_eff.txt',sep='\t') # output from TnP tool
df_relIso_bin_D = pn.read_csv('RelIso_bin_D.txt',sep='\t') # calculated relIso bin centers, average relIso distribution value in the particular (eta,pT,relIso) bin for DATA. Column format -> ||| eta bin (1,2,3) | pT bin (7,10,15) | iso bin left boundary (value) | average relIso value in the particular pT/eta/iso bin (value) |||
df_relIso_bin_M = pn.read_csv('RelIso_bin_M.txt',sep='\t') # calculated relIso bin centers, average relIso distribution value in the particular (eta,pT,relIso) bin for MC.

z_avg_relIso_bins = [0.49,0.36,0.26,0.26,0.20,0.17,0.14,0.12,0.10] #Inner barrel 3 pT bins, outer barrel 3 pT bins, endcap 3 pT bins

#df_eff["d_eff_err_full"] = ( (df_eff["d_eff_n_err"])**2 + (df_eff["d_eff_b"]-df_eff["d_eff_n"])**2 + (df_eff["d_eff_s"]-df_eff["d_eff_n"])**2 + (df_eff["d_eff_t"]-df_eff["d_eff_n"])**2 )**0.5
#df_eff["m_eff_err_full"] = ( (df_eff["m_eff_err"])**2 + (df_eff["m_eff_t"]-df_eff["m_eff"])**2 )**0.5

def fit_func(x,a,b):
    return a*x*x+b
bins = ['5_7_inner','7_10_inner','10_15_inner','5_7_outer','7_10_outer','10_15_outer','5_7_endcap','7_10_endcap','10_15_endcap']
eta_bin = ['0.0-0.8 $|\eta|$','0.0-0.8 $|\eta|$','0.0-0.8 $|\eta|$','0.8-1.44 $|\eta|$','0.8-1.44 $|\eta|$','0.8-1.44 $|\eta|$','1.44-2.5 $|\eta|$','1.44-2.5 $|\eta|$','1.44-2.5 $|\eta|$']
pT_bin = ['5-7 GeV','7-10 GeV','10-15 GeV','5-7 GeV','7-10 GeV','10-15 GeV','5-7 GeV','7-10 GeV','10-15 GeV']
# inner 5-7 GeV
xdata = np.linspace(0, 1, 100)

for i in range(9): # 3 eta bins x 3pT bins == 9 bins in total

    if i < 3 : # Inner Barrel pT bins, relIso bin center value
        x_coord_MC = df_relIso_bin_M['relISO'][(11*i+0):(11*i+4)] # if we want to plot all 9 plots, have to reformat relIso columns
        x_coord_DATA = df_relIso_bin_D['relISO'][(11*i+0):(11*i+4)]

    if i >= 3 or i < 6  : # Outer Barrel pT bins, relIso bin center value
        x_coord_MC = df_relIso_bin_M['relISO'][(11*(i-3)+4):(11*(i-3)+8)] # if we want to plot all 9 plots, have to reformat relIso columns
        x_coord_DATA = df_relIso_bin_D['relISO'][(11*(i-3)+4):(11*(i-3)+8)]

    if i>=6 : # Endcap pT bins, relIso bin center value
        x_coord_MC = df_relIso_bin_M['relISO'][(11*(i-6)+8):(11*(i-6)+11)] # if we want to plot all 9 plots, have to reformat relIso columns
        x_coord_DATA = df_relIso_bin_D['relISO'][(11*(i-6)+8):(11*(i-6)+11)]
    
    # Endcap has 3 relIso bins instead of 4. 
    if i == 6:
        y_coord_MC = df_eff['m_eff'][(4*i+0):(4*i+3)]
        y_coord_MC_Err = df_eff["m_eff_err"][(4*i+0):(4*i+3)]

        y_coord_DATA = df_eff['d_eff_mean'][(4*i+0):(4*i+3)]
        y_coord_DATA_Err=df_eff['d_eff_mean_err'][(4*i+0):(4*i+3)]

    elif  i==7:
        y_coord_MC = df_eff['m_eff'][(4*i-1):(4*i+2)]
        y_coord_MC_Err = df_eff["m_eff_err"][(4*i-1):(4*i+2)]

        y_coord_DATA = df_eff['d_eff_mean'][(4*i-1):(4*i+2)]
        y_coord_DATA_Err=df_eff['d_eff_mean_err'][(4*i-1):(4*i+2)]

    elif  i==8:
        y_coord_MC = df_eff['m_eff'][(4*i-2):(4*i+1)]
        y_coord_MC_Err = df_eff["m_eff_err"][(4*i-2):(4*i+1)]

        y_coord_DATA = df_eff['d_eff_mean'][(4*i-2):(4*i+1)]
        y_coord_DATA_Err=df_eff['d_eff_mean_err'][(4*i-2):(4*i+1)]
    
    else:
        y_coord_MC = df_eff['m_eff'][(4*i+0):(4*i+4)]
        y_coord_MC_Err = df_eff["m_eff_err"][(4*i+0):(4*i+4)]

        y_coord_DATA = df_eff['d_eff_mean'][(4*i+0):(4*i+4)]
        y_coord_DATA_Err=df_eff['d_eff_mean_err'][(4*i+0):(4*i+4)]


    popt_MC, pcov_MC = curve_fit(fit_func, x_coord_MC, y_coord_MC, sigma=y_coord_MC_Err,absolute_sigma=True)
    popt_DATA, pcov_DATA = curve_fit(fit_func, x_coord_DATA, y_coord_DATA, sigma=y_coord_DATA_Err,absolute_sigma=True)

    d_eff = (z_avg_relIso_bins[i]**2*popt_DATA[0]+popt_DATA[1])
    d_eff_err =(z_avg_relIso_bins[i]**4*pcov_DATA[0][0]+pcov_DATA[1][1]+2*z_avg_relIso_bins[i]**2*pcov_DATA[0][1])**0.5
    mc_eff = (z_avg_relIso_bins[i]**2*popt_MC[0]+popt_MC[1])
    mc_eff_err = (z_avg_relIso_bins[i]**4*pcov_MC[0][0]+pcov_MC[1][1]+2*z_avg_relIso_bins[i]**2*pcov_MC[0][1])**0.5
    Sf = d_eff/mc_eff
    Sf_err=Sf*((d_eff_err/d_eff)**2+(mc_eff_err/mc_eff)**2)**0.5

    plt.errorbar(x_coord_MC, y_coord_MC,yerr=y_coord_MC_Err,color='g', label='MC',fmt='o',capsize=5,ms=5,mew=2,alpha=0.8)
    plt.errorbar(x_coord_DATA, y_coord_DATA,yerr=y_coord_DATA_Err,color='r', label='DATA',fmt='o',capsize=5,ms=5,mew=2,alpha=0.8)
    plt.plot(xdata, fit_func(xdata, *popt_MC), 'g--',label='Fit MC')
    plt.plot(xdata, fit_func(xdata, *popt_DATA), 'r--',label='Fit DATA')
    plt.axvline(z_avg_relIso_bins[i],ymin=0,ymax=1.05,ls='--',color='gray',lw=2,alpha=0.8, label='Z relIso {0}'.format(z_avg_relIso_bins[i]),visible=True)
    plt.xlabel('Relative Isolation')
    plt.ylabel('Efficiency')
    plt.text(0.00,1.06,'$\\bf{CMS}$',fontsize = 20)
    plt.text(0.10,1.06,'Preliminary',fontsize = 20)
    plt.title('67.9 $fb^{-1}$ (13 TeV) 2018', loc = 'right', fontsize= 20)
    plt.ylim(0,1.05)
    plt.grid()
    plt.text(0.01,1.0,r'$SF: {0:.3f}\pm{1:.3f}$'.format(Sf,Sf_err),fontweight='bold',fontsize=20) # SF
    plt.text(0.01,0.96,'$\epsilon_{d}$',fontsize=20) # data eff
    plt.text(0.01,0.92,'$\epsilon_{m}$',fontsize=20) # mc eff
    plt.text(0.06,0.96,r'$: {0:.3f}\pm{1:.3f}$'.format(d_eff,d_eff_err),fontsize=20) # data eff
    plt.text(0.06,0.92,r'$: {0:.3f}\pm{1:.3f}$'.format(mc_eff,mc_eff_err),fontsize=20) # mc eff
    plt.text(0.06,0.1,eta_bin[i],fontsize=20) # 
    plt.text(0.06,0.05,pT_bin[i],fontsize=20) #
    plt.xlim(-0.01,1.0)
    plt.legend()
    plt.savefig('/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/DPS_plots/eff_extrapol_{0}_RMS.png'.format(bins[i]))
    plt.savefig('/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/DPS_plots/eff_extrapol_{0}_RMS.pdf'.format(bins[i]))
    plt.show()

