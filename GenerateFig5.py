#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril
    
This script generates a figure in the manuscript from raw experimental data.

Figure 5 Antibiotic activity of benazepril against M. smegmatis
"""

#%% Import packages

# Data packages
import numpy as np
import pandas as pd
import json

# Graphics packages
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 5/'
Fig_meta['Style path'] = 'Style files/'

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)

#%% Read data files used to make the plots

# Load the Msmeg growth curves
Msmeg_IC50_DF = pd.read_csv(Fig_meta['Data path'] + 'Msmeg IC50 Curves.csv', header=[0], index_col=[0, 1, 2])

# Load the Msmeg Dala rescue data
Msmeg_Dala_DF = pd.read_csv(Fig_meta['Data path']  + 'Msmeg Dala Rescue.csv', header=[0], index_col=[0, 1])

#%% Lay out the figure

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] +  'Style Fig 5.mplstyle')

# Arrange the panels using GridSpec
Ax_handle = [None] * 2
Fig_handle = plt.figure(tight_layout=False)

Fig_meta['Rows'] = 3
Fig_meta['Cols'] = 5 
Fig_meta['Height ratios'] = [0.5,2,0.7]
Fig_meta['Width ratios'] = [1,4,1,2,0.5]

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle, 
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'])
                                     
Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1, 1]) # Panel A
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[1, 3]) # Panel B


# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]
    
    Fig_handle.text(Text_X, Text_Y, Label, weight='bold', fontsize=8,
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

add_panel_label(Ax_handle[0], 'A', (-0.1, 0.1))
add_panel_label(Ax_handle[1], 'B', (-0.1, 0.1))

#%% Collect the Msmeg IC50 growth data

Panel_Msmeg = {}

Panel_Msmeg['DMSO'] = {}
Panel_Msmeg['DMSO']['DataFrame'] = Msmeg_IC50_DF.loc['DMSO (%)']
Panel_Msmeg['DMSO']['Concentrations'] = np.sort(Panel_Msmeg['DMSO']['DataFrame'].index.remove_unused_levels().levels[1].values.astype(float))
Panel_Msmeg['DMSO']['OD'] = Msmeg_IC50_DF.loc['DMSO (%)'].values.astype(float)
Panel_Msmeg['DMSO']['Mean'] = np.mean(Panel_Msmeg['DMSO']['OD'], axis=1)
Panel_Msmeg['DMSO']['Std'] = np.std(Panel_Msmeg['DMSO']['OD'], axis=1)
Panel_Msmeg['DMSO']['95 conf'] = (1.96/np.sqrt(12)) * np.std(Panel_Msmeg['DMSO']['OD'], axis=1)

Panel_Msmeg['Benazeprilat'] = {}
Panel_Msmeg['Benazeprilat']['DataFrame'] = Msmeg_IC50_DF.loc['Benazeprilat (mM)']
Panel_Msmeg['Benazeprilat']['Concentrations'] = np.sort(Panel_Msmeg['Benazeprilat']['DataFrame'].index.remove_unused_levels().levels[1].values.astype(float))
Panel_Msmeg['Benazeprilat']['OD'] = Msmeg_IC50_DF.loc['Benazeprilat (mM)'].values.astype(float)
Panel_Msmeg['Benazeprilat']['Mean'] = np.mean(Panel_Msmeg['Benazeprilat']['OD'], axis=1)
Panel_Msmeg['Benazeprilat']['Std'] = np.std(Panel_Msmeg['Benazeprilat']['OD'], axis=1)
Panel_Msmeg['Benazeprilat']['95 conf'] = (1.96/np.sqrt(12)) * np.std(Panel_Msmeg['Benazeprilat']['OD'], axis=1)

Panel_Msmeg['Benazepril'] = {}
Panel_Msmeg['Benazepril']['DataFrame'] = Msmeg_IC50_DF.loc['Benazepril (mM)']
Panel_Msmeg['Benazepril']['Concentrations'] = np.sort(Panel_Msmeg['Benazepril']['DataFrame'].index.remove_unused_levels().levels[1].values.astype(float))
Panel_Msmeg['Benazepril']['OD'] = Msmeg_IC50_DF.loc['Benazepril (mM)'].values.astype(float)
Panel_Msmeg['Benazepril']['Mean'] = np.mean(Panel_Msmeg['Benazepril']['OD'], axis=1)
Panel_Msmeg['Benazepril']['Std'] = np.std(Panel_Msmeg['Benazepril']['OD'], axis=1)
Panel_Msmeg['Benazepril']['95 conf'] = (1.96/np.sqrt(12)) * np.std(Panel_Msmeg['Benazepril']['OD'], axis=1)

Panel_Msmeg['DCS'] = {}
Panel_Msmeg['DCS']['DataFrame'] = Msmeg_IC50_DF.loc['DCS (mM)']
Panel_Msmeg['DCS']['Concentrations'] = np.sort(Panel_Msmeg['DCS']['DataFrame'].index.remove_unused_levels().levels[1].values.astype(float))
Panel_Msmeg['DCS']['OD'] = Msmeg_IC50_DF.loc['DCS (mM)'].values.astype(float)
Panel_Msmeg['DCS']['Mean'] = np.mean(Panel_Msmeg['DCS']['OD'], axis=1)
Panel_Msmeg['DCS']['Std'] = np.std(Panel_Msmeg['DCS']['OD'], axis=1)
Panel_Msmeg['DCS']['95 conf'] = (1.96/np.sqrt(12)) * np.std(Panel_Msmeg['DCS']['OD'], axis=1)

#%% Plot the Msmeg IC50 growth data

Ax = Ax_handle[0]

# Style the lines
Panel_Msmeg['DMSO']['Line properties'] = {'color':KellyColors['Buff']}
Panel_Msmeg['Benazeprilat']['Line properties'] = {'color':KellyColors['Light pink']}
Panel_Msmeg['Benazepril']['Line properties'] = {'color':KellyColors['Dark pink']}
Panel_Msmeg['DCS']['Line properties'] = {'color':KellyColors['Medium blue']}

Ax.errorbar(Panel_Msmeg['DMSO']['Concentrations'], Panel_Msmeg['DMSO']['Mean'], yerr = Panel_Msmeg['DMSO']['95 conf'], **Panel_Msmeg['DMSO']['Line properties'])
Ax.errorbar(Panel_Msmeg['Benazeprilat']['Concentrations'], Panel_Msmeg['Benazeprilat']['Mean'], yerr = Panel_Msmeg['Benazeprilat']['95 conf'], **Panel_Msmeg['Benazeprilat']['Line properties'] )
Ax.errorbar(Panel_Msmeg['Benazepril']['Concentrations'], Panel_Msmeg['Benazepril']['Mean'], yerr = Panel_Msmeg['Benazepril']['95 conf'], **Panel_Msmeg['Benazepril']['Line properties'])
Ax.errorbar(Panel_Msmeg['DCS']['Concentrations'], Panel_Msmeg['DCS']['Mean'], yerr = Panel_Msmeg['DCS']['95 conf'], **Panel_Msmeg['DCS']['Line properties'])

# Label and style the axes
Ax.set_xlabel('Concentration (mM or %)')
Ax.set_xscale('Log')
Ax.set_ylim([0,0.45])
Ax.set_yticks([0,0.1,0.2,0.3,0.4])
Ax.set_ylabel('Growth (OD)')

Ax.legend(['DMSO (%)', 'Benazeprilat', 'Benazepril', 'DCS'], frameon=False, loc='lower left', handlelength=1, labelspacing=0.2)

#%% Collect the Msmeg Dala rescue data

Panel_Rescue = {}
Panel_Rescue['Mean'] = np.mean(Msmeg_Dala_DF.values, axis=1)
Panel_Rescue['Std'] = np.std(Msmeg_Dala_DF.values, axis=1)

#%% Plot the Msmeg Dala rescue data

Ax = Ax_handle[1]

Panel_Rescue['Positions'] = [1, 2, 4, 5, 7, 8]

Panel_Rescue['Bars'] = Ax.bar(Panel_Rescue['Positions'], Panel_Rescue['Mean'], yerr=Panel_Rescue['Std'])

# Style the Dala rescue bar graph
for idx in [0,1]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Buff'])
for idx in [2,3]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Medium blue'])
for idx in [4,5]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Dark pink'])
for idx in [0,2,4]: Panel_Rescue['Bars'][idx].set_hatch('////')

# Label and style the axes
Ax.set_xticks([1.5, 4.5, 7.5])
Ax.set_xticklabels(['No\nDrug', 'DCS', 'Benaz.'])

Ax.set_ylim([0,0.45])
Ax.set_yticks([0,0.1,0.2,0.3,0.4])
Ax.set_ylabel('Growth (OD)')

# Make a custom legend
Panel_Rescue['Hatch patch'] = mpatches.Patch(facecolor='#FFFFFF', hatch='////', label='+ D-ala')
Ax.legend(handles = [Panel_Rescue['Hatch patch']], bbox_to_anchor=(0.5, 0.95), loc='center', handlelength=1.55, frameon=False)

#%% Save the figure as a pdf

plt.savefig('Fig 5.pdf')



