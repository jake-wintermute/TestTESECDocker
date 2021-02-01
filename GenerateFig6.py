#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril
    
This script generates a figure in the manuscript from raw experimental data.

Figure 6 Benazepril inhibits Mtb ALR in vitro. 
"""

#%% Import packages

# Data packages
import numpy as np
import pandas as pd
import json
from scipy.optimize import curve_fit

# Graphics packages
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 6/'
Fig_meta['Style path'] = 'Style files/'

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)
Meta = {};

#%% Read data files used to make the plots

HW_DF  = pd.read_csv(Fig_meta['Data path'] + 'In Vitro Hanes Woolf.csv', index_col=[0,1], header=[0,1])
MM_DF  = pd.read_csv(Fig_meta['Data path'] + 'In Vitro Michaelis Menten.csv', index_col=[0,1], header=[0,1])

#%% Lay out the figure

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] +  'Style Fig 6.mplstyle')

# Arrange the panels using GridSpec
Fig_meta['Rows'] = 5
Fig_meta['Cols'] = 13 
Fig_meta['Height ratios'] = [1,1,1.5,1,1]
Fig_meta['Width ratios'] = [1,1,1,1,1,1,1,1,1,1,1,1,0.5]

Ax_handle = [None] * 6
Fig_handle = plt.figure(tight_layout=False)

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle, 
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'])
                                     
Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1, 1:3]) # Panel A top
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[3, 1:3]) # Panel A bottom
Ax_handle[2] = Fig_handle.add_subplot(GS_handle[1:4, 4:8]) # Panel B
Ax_handle[3] = Fig_handle.add_subplot(GS_handle[1, 10:12]) # Panel C top
Ax_handle[4] = Fig_handle.add_subplot(GS_handle[3, 10:12]) # Panel C bottom

# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]
    
    Fig_handle.text(Text_X, Text_Y, Label, weight='bold', fontsize=8,
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

add_panel_label(Ax_handle[0], 'A', (-0.05, 0.08))
add_panel_label(Ax_handle[2], 'B', (-0.05, 0.08))
add_panel_label(Ax_handle[3], 'C', (-0.05, 0.08))

#%% Collect the Hanes-Woolf into a dict for numerical analysis

HW_Dict = {}

HW_Dict['Benazepril (mM)'] = np.sort(HW_DF.columns.remove_unused_levels().levels[1].values.astype(float))
HW_Dict['Alanine (mM)'] = np.sort(HW_DF.index.remove_unused_levels().levels[0].values.astype(float))
HW_Dict['Replicates'] =  np.sort(HW_DF.index.remove_unused_levels().levels[1].values.astype(str))

#%% Build the Hanes-Woolf data into a multidimensional array

HW_Dict['Vf matrix'] = np.ndarray([len(HW_Dict['Benazepril (mM)']), len(HW_Dict['Alanine (mM)']), len(HW_Dict['Replicates'])])
for i, Benaz_mm in enumerate(HW_Dict['Benazepril (mM)']):
    for j, Ala_mm in enumerate(HW_Dict['Alanine (mM)']):
        for k, Rep in enumerate(HW_Dict['Replicates']):
            
            HW_Dict['Vf matrix'][i,j,k] =  HW_DF.loc[(Ala_mm, Rep), ('Benazepril (mM)', str(Benaz_mm))]   

HW_Dict['Vf mean'] = np.mean(HW_Dict['Vf matrix'], axis=2)
HW_Dict['Vf std'] = np.std(HW_Dict['Vf matrix'], axis=2)
HW_Dict['Vf 95 conf'] = 1.97 * (HW_Dict['Vf std'] / np.sqrt(3))

# Clean up the workspace
del i, j, k, Benaz_mm, Ala_mm

#%% Fit parameter values for the Hanes-Woolf data

Fit_MM = {};
Fit_MM['Vf matrix']= HW_Dict['Vf matrix'] 

Fit_MM['Benazepril (mM)'] = HW_Dict['Benazepril (mM)'] 
Fit_MM['Alanine (mM)'] = HW_Dict['Alanine (mM)']

# Replicate the alanine concentrations for all replicates
Fit_MM['Alanine reps'] = np.tile(Fit_MM['Alanine (mM)'], (3,1)).T

# Perform the curve fitting to a Michaelis-Menten model
Cfit_Params = {}
Cfit_Params['bounds'] = ([0, 0], [200, 20])

def Michaels_Menten_func(S, Vmax, Km):
    return (Vmax * S) / (Km + S)

# Fit the parameters for each benazepril concentration
Fit_MM['Fit params'] = {}
Fit_MM['Fit params']['Solver out'] = [None] * len(Fit_MM['Benazepril (mM)'])
Fit_MM['Fit params']['Solver cov'] = [None] * len(Fit_MM['Benazepril (mM)'])

for i in range(len(Fit_MM['Benazepril (mM)'])):
    Fit_MM['Fit params']['Solver out'][i], Fit_MM['Fit params']['Solver cov'][i] = curve_fit(Michaels_Menten_func,
                                                    np.ndarray.flatten(Fit_MM['Alanine reps']),
                                                    np.ndarray.flatten(Fit_MM['Vf matrix'][i, :, :]),
                                                    **Cfit_Params)

Fit_MM['Fit params']['Vmax'] =  np.array([Solver_Output[0] for Solver_Output in Fit_MM['Fit params']['Solver out']])
Fit_MM['Fit params']['Km'] =  np.array([Solver_Output[1] for Solver_Output in Fit_MM['Fit params']['Solver out']])

# Produce the best-fit curve lines
Fit_MM['Fit Vf curves'] = {}
Fit_MM['Fit Vf curves']['Alanine (mM)'] = np.arange(0.1, max(Fit_MM['Alanine (mM)']), 0.1)
Fit_MM['Fit Vf curves']['Vf'] = [None] * len(Fit_MM['Benazepril (mM)'])

for b in range(len(Fit_MM['Benazepril (mM)'])):
    Fit_MM['Fit Vf curves']['Vf'][b] = Michaels_Menten_func(Fit_MM['Fit Vf curves']['Alanine (mM)'], Fit_MM['Fit params']['Vmax'][b], Fit_MM['Fit params']['Km'][b])

# Extract the parameter confidence intevals from the coviariance matrix
Fit_MM['Fit params']['Vmax var'] = [CovMat[0,0] for CovMat in Fit_MM['Fit params']['Solver cov']]
Fit_MM['Fit params']['Km var']  = [CovMat[1,1] for CovMat in Fit_MM['Fit params']['Solver cov']]

Fit_MM['Fit params']['Vmax 95 conf'] = (1.97/np.sqrt(3)) * np.sqrt(Fit_MM['Fit params']['Vmax var'])
Fit_MM['Fit params']['Km 95 conf']  = (1.97/np.sqrt(3)) * np.sqrt(Fit_MM['Fit params']['Km var'])      

# Clean up the workspace
del i, b

#%% Collect the data for the Michelis-Menten curves

Panel_MM = {}

Panel_MM['Benazepril'] = {}
Panel_MM['DCS'] = {}

Panel_MM['Benazepril']['Concentration (mM)'] = np.sort(MM_DF.loc['Benazepril'].index.values.astype(float))
Panel_MM['DCS']['Concentration (mM)'] = np.sort(MM_DF.loc['DCS'].index.values.astype(float))

Panel_MM['Benazepril']['Vf all'] = MM_DF.loc['Benazepril'].values
Panel_MM['Benazepril']['Vf mean'] = np.mean(Panel_MM['Benazepril']['Vf all'], axis=1)
Panel_MM['Benazepril']['Vf 95 conf'] = (1.97 / np.sqrt(3)) * np.std(Panel_MM['Benazepril']['Vf all'], axis=1)

Panel_MM['DCS']['Vf all'] = MM_DF.loc['DCS'].values
Panel_MM['DCS']['Vf mean'] = np.mean(Panel_MM['DCS']['Vf all'], axis=1)
Panel_MM['DCS']['Vf 95 conf'] = (1.97 / np.sqrt(3)) * np.std(Panel_MM['DCS']['Vf all'], axis=1)

#%% Plot the data for the Michaelis-Menten curves

Ax = Ax_handle[0]
Panel_MM['Benazepril']['Line properties'] = {'color':KellyColors['Dark pink']}

Ax.errorbar(Panel_MM['Benazepril']['Concentration (mM)'], Panel_MM['Benazepril']['Vf mean'], yerr = Panel_MM['Benazepril']['Vf 95 conf'], **Panel_MM['Benazepril']['Line properties'] )

# Label and style the axes
Ax.set_yticks([0,50,100])
Ax.set_xscale('log')
Ax.set_xlabel('Benazepril (mM)')
Ax.set_ylabel('ALR V$_F$')


Ax = Ax_handle[1]
Panel_MM['DCS']['Line properties'] = {'color':KellyColors['Medium blue']}

Ax.errorbar(Panel_MM['DCS']['Concentration (mM)'], Panel_MM['DCS']['Vf mean'], yerr = Panel_MM['DCS']['Vf 95 conf'], **Panel_MM['DCS']['Line properties'] )

# Label and style the axes
Ax.set_yticks([0,50,100])
Ax.set_xscale('log')
Ax.set_xlim([1E-4, 1])
Ax.set_xlabel('DCS (mM)')
Ax.set_ylabel('ALR V$_F$')

#%% Collect the data for the Michaelis-Menten series

Panel_MM_Series = {};
Panel_MM_Series['Benazepril (mM)'] = Fit_MM['Benazepril (mM)']
Panel_MM_Series['Alanine (mM)'] = Fit_MM['Alanine (mM)']
Panel_MM_Series['Fit Vf curves'] = Fit_MM['Fit Vf curves']
Panel_MM_Series['Vf matrix'] = Fit_MM['Vf matrix']

#%% Plot the Michaelis-Menten series

Ax = Ax_handle[2]

# Define the color and style of the lines
Panel_MM_Series['Colormap'] = plt.cm.inferno_r
Panel_MM_Series['Num lines'] = 9
Panel_MM_Series['Line colors'] = Panel_MM_Series['Colormap'](np.linspace(0, 1, Panel_MM_Series['Num lines']))

Panel_MM_Series['Marker properties'] = {'alpha':0.5, 'markersize':4, 'markeredgecolor':'none', 'marker':'o', 'linestyle':'none'}

# Plot the best-fit Michaelis-Menten curves
for b in range(len(Panel_MM_Series['Benazepril (mM)'])):
    Ax.plot(Panel_MM_Series['Fit Vf curves']['Alanine (mM)'], Panel_MM_Series['Fit Vf curves']['Vf'][b], color=Panel_MM_Series['Line colors'][b])
    
# Plot the individual points of the kinetic curves
for b in range(len(Panel_MM_Series['Benazepril (mM)'])):
        for r in range(3):
            Ax.plot(Panel_MM_Series['Alanine (mM)'], Panel_MM_Series['Vf matrix'][b, :, r], color=Panel_MM_Series['Line colors'][b], **Panel_MM_Series['Marker properties'])

# Label and style the axes
Ax.set_xlabel('Alanine (mM)')
Ax.set_ylabel('ALR V$_F$')

# Make an axis for the colorbar
Panel_MM_Series['Color norm']  = mpl.colors.Normalize(vmin=0, vmax=1)

Panel_MM_Series['Colorbar ax'] = inset_axes(Ax,
                    width="100%",
                    height="100%",
                    bbox_to_anchor=(1.08, 0.1, 0.05, 0.8),
                    bbox_transform=Ax.transAxes,
                    borderpad=0)

Panel_MM_Series['Colorbar handle'] = mpl.colorbar.ColorbarBase(
    Panel_MM_Series['Colorbar ax'],
    cmap=Panel_MM_Series['Colormap'],
    norm=Panel_MM_Series['Color norm'],
    orientation='vertical')

Panel_MM_Series['Colorbar handle'].ax.yaxis.set_label_position('left')
Panel_MM_Series['Colorbar handle'].ax.set_ylabel('Benazepril (mM)')
Panel_MM_Series['Colorbar handle'].ax.invert_yaxis()

#%% Collect the best-fit kinetic parameters

Panel_Params = {};

#Truncate to only the benazepril concetrations that allowed parameter fitting
Panel_Params['Good fits'] = slice(0,-2)

Panel_Params['Benazepril (mM)'] = Fit_MM['Benazepril (mM)'][Panel_Params['Good fits']]
Panel_Params['Vmax'] = Fit_MM['Fit params']['Vmax'][Panel_Params['Good fits']]
Panel_Params['Km'] = Fit_MM['Fit params']['Km'][Panel_Params['Good fits']]
Panel_Params['Vmax 95 conf'] = Fit_MM['Fit params']['Vmax 95 conf'][Panel_Params['Good fits']]
Panel_Params['Km 95 conf'] = Fit_MM['Fit params']['Km 95 conf'][Panel_Params['Good fits']]

#%% Plot the best-fit kinetic parameters

Ax = Ax_handle[3]

Panel_Params['Km properties'] = {'color':KellyColors['Medium green']}

Ax.errorbar(Panel_Params['Benazepril (mM)'], Panel_Params['Vmax'], yerr=Panel_Params['Vmax 95 conf'], **Panel_Params['Km properties'])

# Label and style the axes
Ax.set_ylim([0,125])
Ax.set_yticks([0,50,100])
Ax.set_xlabel('Benazepril (mM)')
Ax.set_ylabel('ALR V$_{MAX}$')

Ax = Ax_handle[4]

Panel_Params['Vmax properties'] = {'color':KellyColors['Violet']}

Ax.errorbar(Panel_Params['Benazepril (mM)'], Panel_Params['Km'], yerr=Panel_Params['Km 95 conf'], **Panel_Params['Vmax properties'])

# Label and style the axes
Ax.set_ylim([0,10])
Ax.set_yticks([0,5,10])
Ax.set_xlabel('Benazepril (mM)')
Ax.set_ylabel('ALR K$_M$')

#%% Save the figure as a pdf

plt.savefig('Fig 6.pdf')