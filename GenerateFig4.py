#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril

This script generates a figure in the manuscript from raw experimental data.

Figure 4 Chemical-genetic profiles indicate distinct mechanisms of action
"""

#%% Import packages

# Data packages
import numpy as np
import pandas as pd
import json

# Graphics packages
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as mpatches

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 4/'
Fig_meta['Style path'] = 'Style files/'
Fig_meta['Drug names'] = ['DCS', 'Benazepril', 'Amlexanox']

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)

#%% Read data files used to make the plots

# Load the drug screen data
ChemGen_DF = pd.read_csv(Fig_meta['Data path']  + 'TESEC ChemGen Profiles.csv', header=[0], index_col=[0, 1, 2])
Rescue_DF = pd.read_csv(Fig_meta['Data path']  + 'TESEC Dala Rescue.csv', header=[0], index_col=[0, 1])
Efflux_DF = pd.read_csv(Fig_meta['Data path'] + 'TESEC Efflux Effect.csv', header=[0], index_col=[0,1,2])

#%% Build the figure layout

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] +  'Style fig 4.mplstyle')

# Arrange the panels using GridSpec
Ax_handle = [None] * 8
Fig_handle = plt.figure(tight_layout=False)

Fig_meta['Rows'] = 5
Fig_meta['Cols'] = 16
Fig_meta['Height ratios'] = [0.5,1,0.8,1,0.5]
Fig_meta['Width ratios'] = [1.2,1,1,1,1.2,1,1,1,1.2,1,1,1,1.2,1,1,0.5]

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle,
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'])

Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1, 1:4]) # Panel A
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[1, 5:8]) # Panel B
Ax_handle[2] = Fig_handle.add_subplot(GS_handle[1, 9:12]) # Panel C
Ax_handle[3] = Fig_handle.add_subplot(GS_handle[1, 13:15]) # Panel D
Ax_handle[4] = Fig_handle.add_subplot(GS_handle[3, 1:4]) # Panel E
Ax_handle[5] = Fig_handle.add_subplot(GS_handle[3, 5:8]) # Panel F
Ax_handle[6] = Fig_handle.add_subplot(GS_handle[3, 9:12]) # Panel G
Ax_handle[7] = Fig_handle.add_subplot(GS_handle[3, 13:15]) # Panel H

# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]

    Fig_handle.text(Text_X, Text_Y, Label, weight='bold', fontsize=8,
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

add_panel_label(Ax_handle[0], 'A', (-0.05, 0.08))
add_panel_label(Ax_handle[1], 'B', (-0.05, 0.08))
add_panel_label(Ax_handle[2], 'C', (-0.05, 0.08))
add_panel_label(Ax_handle[3], 'D', (-0.05, 0.08))
add_panel_label(Ax_handle[4], 'E', (-0.05, 0.08))
add_panel_label(Ax_handle[5], 'F', (-0.05, 0.08))
add_panel_label(Ax_handle[6], 'G', (-0.05, 0.08))
add_panel_label(Ax_handle[7], 'H', (-0.05, 0.08))

#%% Collect the data for the heatmaps

Panel_Hmap = {}
for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name] = {}

for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['OD mean'] = ChemGen_DF.loc[Drug_Name, 'Mean']
    Panel_Hmap[Drug_Name]['Arabinose concentrations'] = Panel_Hmap[Drug_Name]['OD mean'].columns.values.astype('float')
    Panel_Hmap[Drug_Name]['Drug concentrations'] = (1/1000) * Panel_Hmap[Drug_Name]['OD mean'].index.values.astype('float')

# Limit the heatmap to only the top 8 rows
Panel_Hmap['Benazepril']['Top rows'] = np.array(range(1,8))
Panel_Hmap['DCS']['Top rows'] = np.array(range(1,8))
Panel_Hmap['Amlexanox']['Top rows'] = np.array(range(1,8))

for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['OD mean'] = Panel_Hmap[Drug_Name]['OD mean'].iloc[Panel_Hmap[Drug_Name]['Top rows'], :]
    Panel_Hmap[Drug_Name]['Drug concentrations'] = Panel_Hmap[Drug_Name]['Drug concentrations'][Panel_Hmap[Drug_Name]['Top rows']]

# Add additional points to define the corners of the heatmap boxes
for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['X'] = np.append(Panel_Hmap[Drug_Name]['Arabinose concentrations'], 2*Panel_Hmap[Drug_Name]['Arabinose concentrations'][-1])
    Panel_Hmap[Drug_Name]['Y'] = np.append(Panel_Hmap[Drug_Name]['Drug concentrations'], 2*Panel_Hmap[Drug_Name]['Drug concentrations'][-1])

#%% Style the heatmaps

Panel_Hmap['DCS']['Ax handle'] = Ax_handle[0]
Panel_Hmap['Benazepril']['Ax handle'] = Ax_handle[1]
Panel_Hmap['Amlexanox']['Ax handle'] = Ax_handle[2]

Panel_Hmap['Hmap kwargs'] = {'edgecolors':'w', 'linewidths':0.1,
                                 'cmap':'viridis', 'vmin':0.0, 'vmax':0.4}

Panel_Hmap['X ticks'] = [1E2, 1E5, 1E8]
Panel_Hmap['X label'] = 'Arabinose (nM)'

Panel_Hmap['DCS']['Y label'] = 'DCS (mM)'
Panel_Hmap['Benazepril']['Y label'] = 'Benazepril (mM)'
Panel_Hmap['Amlexanox']['Y label'] = 'Amlexanox (mM)'

# Style the heatmap colorbars
Panel_Hmap['Colorbar kwargs']  = {'orientation':'horizontal'}

# Define the colorbar axes for the heat maps
for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['Colorbar cax'] = make_axes_locatable(Panel_Hmap[Drug_Name]['Ax handle']).append_axes(
        'top', size='15%', pad='30%')

#%% Plot the heatmaps

for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['Pcolor handle'] = Panel_Hmap[Drug_Name]['Ax handle'].pcolor(Panel_Hmap[Drug_Name]['X'] , Panel_Hmap[Drug_Name]['Y'],
                                                                   Panel_Hmap[Drug_Name]['OD mean'], **Panel_Hmap['Hmap kwargs'])

# Label the heatmaps
for Drug_Name in Fig_meta['Drug names']:
    Panel_Hmap[Drug_Name]['Ax handle'].set_xscale('log')
    Panel_Hmap[Drug_Name]['Ax handle'].set_yscale('log')
    Panel_Hmap[Drug_Name]['Ax handle'].set_xticks(Panel_Hmap['X ticks']);
    Panel_Hmap[Drug_Name]['Ax handle'].set_xlabel(Panel_Hmap['X label']);
    Panel_Hmap[Drug_Name]['Ax handle'].set_ylabel(Panel_Hmap[Drug_Name]['Y label']);


Panel_Hmap['DCS']['Ax handle'].set_yticks([1E-2, 1E-1])
Panel_Hmap['Benazepril']['Ax handle'].set_yticks([1E-1, 1E0])
Panel_Hmap['Amlexanox']['Ax handle'].set_yticks([1E-1, 1E0])

# Add the heatmap colorbars
for Drug_Name in Fig_meta['Drug names']:
    plt.colorbar(Panel_Hmap[Drug_Name]['Pcolor handle'], cax=Panel_Hmap[Drug_Name]['Colorbar cax'], **Panel_Hmap['Colorbar kwargs'])
    Panel_Hmap[Drug_Name]['Colorbar cax'].xaxis.set_ticks_position('top')
    Panel_Hmap[Drug_Name]['Colorbar cax'].xaxis.set_label_position('bottom')
    Panel_Hmap[Drug_Name]['Colorbar cax'].set_xlabel('Growth (OD)', labelpad=3)
    Panel_Hmap[Drug_Name]['Colorbar cax'].xaxis.set_tick_params(pad=0)

#%% Collect the data for the dose-response curves

Panel_DR = {}
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name] = {}

for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['OD mean'] = ChemGen_DF.loc[Drug_Name, 'Mean']
    Panel_DR[Drug_Name]['OD Stds'] = ChemGen_DF.loc[Drug_Name, 'Std']
    Panel_DR[Drug_Name]['Arabinose concentrations'] = Panel_DR[Drug_Name]['OD mean'].columns.values.astype('float')
    Panel_DR[Drug_Name]['Drug concentrations'] = (1/1000) * Panel_DR[Drug_Name]['OD mean'].index.values.astype('float')

# Collect the dose-response for the no drug treatment
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['No Drug'] = Panel_DR[Drug_Name]['OD mean'].iloc[0, :]

# Limit the plot to only the top 8 rows
Panel_DR['Benazepril']['Top rows'] = np.array(range(1,8))
Panel_DR['DCS']['Top rows'] = np.array(range(1,8))
Panel_DR['Amlexanox']['Top rows'] = np.array(range(1,8))

for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['OD mean'] = Panel_DR[Drug_Name]['OD mean'].iloc[Panel_DR[Drug_Name]['Top rows'], :]
    Panel_DR[Drug_Name]['OD Stds'] = Panel_DR[Drug_Name]['OD Stds'].iloc[Panel_DR[Drug_Name]['Top rows'], :]
    Panel_DR[Drug_Name]['Drug concentrations'] = Panel_DR[Drug_Name]['Drug concentrations'][Panel_DR[Drug_Name]['Top rows']]
    Panel_DR[Drug_Name]['Log drug concentrations'] = np.log10(Panel_DR[Drug_Name]['Drug concentrations'])

#%% Style the dose-response curves

Panel_DR['DCS']['Ax handle'] = Ax_handle[4]
Panel_DR['Benazepril']['Ax handle'] = Ax_handle[5]
Panel_DR['Amlexanox']['Ax handle'] = Ax_handle[6]

Panel_DR['Errorbar kwargs'] = {'linewidth':1}

Panel_DR['X ticks'] = [1E2, 1E5, 1E8]
Panel_DR['X label'] = 'Arabinose (nM)'

Panel_DR['Y limits'] = [0, 0.5]
Panel_DR['Y ticks'] = [0.2, 0.4]

Panel_DR['DCS']['Y label'] = 'Growth (OD)'
Panel_DR['Benazepril']['Y label'] = 'Growth (OD)'
Panel_DR['Amlexanox']['Y label'] = 'Growth (OD)'

# Set the colors of the lines
Panel_DR['Colormap'] = plt.cm.inferno_r
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Num lines'] = len(Panel_DR[Drug_Name]['Top rows'])
    Panel_DR[Drug_Name]['Line colors'] = Panel_DR['Colormap'](np.linspace(0, 1, Panel_DR[Drug_Name]['Num lines']))
    Panel_DR[Drug_Name]['Ax handle'].set_prop_cycle(color=Panel_DR[Drug_Name]['Line colors'])

# Style the heatmap colorbars
Panel_DR['Colorbar kwargs']  = {'orientation':'horizontal'}

for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Color norm'] = mpl.colors.Normalize(vmin=Panel_DR[Drug_Name]['Log drug concentrations'][0],
                                                             vmax=Panel_DR[Drug_Name]['Log drug concentrations'][-1])

# Define the colorbar axes for the dose-response curves
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Colorbar cax'] = make_axes_locatable(Panel_DR[Drug_Name]['Ax handle']).append_axes(
        'top', size='15%', pad='30%')

# Add the colorbars to the figure
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Colorbar handle'] = mpl.colorbar.ColorbarBase(
        Panel_DR[Drug_Name]['Colorbar cax'],
        cmap=Panel_DR['Colormap'],
        norm=Panel_DR[Drug_Name]['Color norm'],
        orientation='horizontal')

# Add the labels to the dose-response colorbars
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Colorbar cax'].xaxis.set_ticks_position('top')
    Panel_DR[Drug_Name]['Colorbar cax'].xaxis.set_label_position('bottom')

# Label and style the colorbar axes
Panel_DR['DCS']['Colorbar cax'].set_xlabel('DCS (mM)', labelpad=3)
Panel_DR['Benazepril']['Colorbar cax'].set_xlabel('Benazepril (mM)', labelpad=3)
Panel_DR['Amlexanox']['Colorbar cax'].set_xlabel('Amlexanox (mM)', labelpad=3)

Panel_DR['DCS']['Colorbar cax'].xaxis.set_ticks([-2, -1])
Panel_DR['DCS']['Colorbar cax'].xaxis.set_ticklabels(['$10^{-2}$', '$10^{-1}$'])
Panel_DR['DCS']['Colorbar cax'].xaxis.set_tick_params(pad=0)

Panel_DR['Benazepril']['Colorbar cax'].xaxis.set_ticks([-1, 0])
Panel_DR['Benazepril']['Colorbar cax'].xaxis.set_ticklabels(['$10^{-1}$', '$10^{0}$'])
Panel_DR['Benazepril']['Colorbar cax'].xaxis.set_tick_params(pad=0)

Panel_DR['Amlexanox']['Colorbar cax'].xaxis.set_ticks([-1, 0])
Panel_DR['Amlexanox']['Colorbar cax'].xaxis.set_ticklabels(['$10^{-1}$', '$10^{0}$'])
Panel_DR['Amlexanox']['Colorbar cax'].xaxis.set_tick_params(pad=0)

#%% Plot the dose-response curves

for Drug_Name in Fig_meta['Drug names']:
    for r in range(len(Panel_DR[Drug_Name]['Top rows'])):
        Panel_DR[Drug_Name]['Ax handle'].errorbar(Panel_DR[Drug_Name]['Arabinose concentrations'],
                                                  Panel_DR[Drug_Name]['OD mean'].iloc[r, :], yerr=Panel_DR[Drug_Name]['OD Stds'].iloc[r, :],
                                                  **Panel_DR['Errorbar kwargs'])

# Plot the no drug treatment
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Ax handle'].plot(Panel_DR[Drug_Name]['Arabinose concentrations'], Panel_DR[Drug_Name]['No Drug'],
                                          color='black', linestyle='dashed')

# Label the dose-response curves
for Drug_Name in Fig_meta['Drug names']:
    Panel_DR[Drug_Name]['Ax handle'].set_xscale('log')
    Panel_DR[Drug_Name]['Ax handle'].set_xticks(Panel_DR['X ticks']);
    Panel_DR[Drug_Name]['Ax handle'].set_xlabel(Panel_DR['X label']);
    Panel_DR[Drug_Name]['Ax handle'].set_yticks(Panel_DR['Y ticks']);
    Panel_DR[Drug_Name]['Ax handle'].set_ylim(Panel_DR['Y limits']);
    Panel_DR[Drug_Name]['Ax handle'].set_ylabel(Panel_DR[Drug_Name]['Y label']);


#%% Collect the data for the Dala growth rescue

Panel_Rescue = {}

Panel_Rescue['OD mean'] = np.mean(Rescue_DF.values, axis=1)
Panel_Rescue['OD std'] = np.std(Rescue_DF.values, axis=1)

#%% Plot the Dala growth rescue

Ax = Ax_handle[3]

Panel_Rescue['Positions'] = [1, 2, 4, 5, 7, 8]

Panel_Rescue['Bars'] = Ax.bar(Panel_Rescue['Positions'], Panel_Rescue['OD mean'], yerr=Panel_Rescue['OD std'])

# Style the Dala rescue bar graph
for idx in [0,1]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Buff'])
for idx in [2,3]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Medium blue'])
for idx in [4,5]: Panel_Rescue['Bars'][idx].set_facecolor(KellyColors['Dark pink'])
for idx in [0,2,4]: Panel_Rescue['Bars'][idx].set_hatch('////')

# Label and style the axes
Ax.set_xticks([1.5, 4.5, 7.5])
Ax.set_xticklabels(['No\nDrug', 'DCS', 'Benaz.'])
Ax.set_ylabel('Growth (OD)')

# Make the custom legend for the hatchmarked bars
Panel_Rescue['Hatch patch'] = mpatches.Patch(facecolor='#FFFFFF', alpha=0, hatch='////', label='+ D-ala')
Ax.legend(handles = [Panel_Rescue['Hatch patch']], bbox_to_anchor=(0.5, 1.05), loc='center', handlelength=1.55, frameon=False)

#%% Collect the data for the efflux effect

Panel_Efflux = {}

# Convert the arabinose levels to a log scale
Panel_Efflux['Arabinose levels'] = Efflux_DF.columns.values.astype(float)
Panel_Efflux['Arabinose levels'][0] = 10;
Panel_Efflux['Arabinose levels'] = np.log10(Panel_Efflux['Arabinose levels'])

Panel_Efflux['All strains'] = ['TESEC Mtb ALR Efflux+', 'TESEC Mtb ALR']
Panel_Efflux['All treatments'] = ['Benazepril', 'Benazeprilat', 'DCS']

for strain in Panel_Efflux['All strains']:
    Panel_Efflux[strain] = {}
    for treatment in Panel_Efflux['All treatments']:
        Panel_Efflux[strain][treatment] = {}
        Panel_Efflux[strain][treatment]['All ODs'] = Efflux_DF.loc[strain].loc[treatment].values.astype(float)
        Panel_Efflux[strain][treatment]['OD mean'] = np.mean( Panel_Efflux[strain][treatment]['All ODs'], axis = 0)
        Panel_Efflux[strain][treatment]['OD std'] = np.std(Panel_Efflux[strain][treatment]['All ODs'], axis = 0)

#%% Plot the efflux effect

Ax = Ax_handle[7]

Panel_Efflux['Colors'] = (KellyColors['Medium blue'], KellyColors['Dark pink'])

strain = 'TESEC Mtb ALR Efflux+'
Ax.errorbar(Panel_Efflux['Arabinose levels'], Panel_Efflux[strain]['DCS']['OD mean'],
            yerr=Panel_Efflux[strain]['DCS']['OD std'], color=Panel_Efflux['Colors'][0], linestyle=':')
Ax.errorbar(Panel_Efflux['Arabinose levels'], Panel_Efflux[strain]['Benazepril']['OD mean'],
            yerr=Panel_Efflux[strain]['Benazepril']['OD std'], color=Panel_Efflux['Colors'][1], linestyle=':')

strain = 'TESEC Mtb ALR'
Ax.errorbar(Panel_Efflux['Arabinose levels'], Panel_Efflux[strain]['DCS']['OD mean'],
            yerr=Panel_Efflux[strain]['DCS']['OD std'], color=Panel_Efflux['Colors'][0])
Ax.errorbar(Panel_Efflux['Arabinose levels'], Panel_Efflux[strain]['Benazepril']['OD mean'],
            yerr=Panel_Efflux[strain]['Benazepril']['OD std'], color=Panel_Efflux['Colors'][1])

# Label and style the axes
Ax.set_xlabel('Arabinose (mM)')
Ax.set_xticks([3, 5, 7])
Ax.set_xticklabels(['$10^{3}$', '$10^{5}$', '$10^{7}$'])

Ax.set_ylabel('Growth (OD)')
Ax.set_ylim(0, 0.5)
Ax.set_yticks([0.2, 0.4])

# Add a custom legend
Panel_Efflux['Legend Lines'] = [mpl.lines.Line2D([0], [0], color=KellyColors['Medium blue'], linestyle='-'),
                                mpl.lines.Line2D([0], [0], color=KellyColors['Dark pink'], linestyle='-'),
                                mpl.lines.Line2D([0], [0], color=KellyColors['Dark gray'], linestyle='-'),
                                mpl.lines.Line2D([0], [0], color=KellyColors['Dark gray'], linestyle=':')]

Panel_Efflux['Legend ax'] = make_axes_locatable(Ax).append_axes('top', size='15%', pad='30%')
Panel_Efflux['Legend ax'].axis('off')

Panel_Efflux['Legend ax'].legend(Panel_Efflux['Legend Lines'], ['DCS', 'Benaz.', 'Δ Efflux', 'Efflux +'],
          loc='center', handlelength=0.8, frameon=False, bbox_to_anchor=(0.5, 0.5),
          ncol=2, columnspacing=0.5)

#%% Save the figure as a pdf

plt.savefig('Fig 4.pdf')
