#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril
    
This script generates a figure in the manuscript from raw experimental data.

Figure 3 Chemical-genetic drug sensitivity is significantly altered by target overexpression. 
"""

#%% Import packages

# Data packages
import numpy as np
from scipy import stats
import pickle
import json

# Graphics packages
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib.colors import ListedColormap

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 3/'
Fig_meta['Style path'] = 'Style files/'

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)
    
#%% Read data files used to make the plots

with open(Fig_meta['Data path'] + 'TESEC ALR Screen Dict 384.pickle', 'rb') as Pickle_File_Handle:
    Screen_Dict_384  = pickle.load(Pickle_File_Handle)

#%% Build the figure layout

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] + 'Style Fig 3.mplstyle')

# Arrange the panels using GridSpec
Fig_handle = plt.figure(tight_layout=False)
Ax_handle = [None] * 7

Fig_meta['Rows'] = 6
Fig_meta['Cols'] = 17 
Fig_meta['Height ratios'] = [0.7,1,1,1,1,0.8]
Fig_meta['Width ratios'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle, 
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'] )

Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1:3, 2:4]) # Panel A
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[1:3, 5:7]) # Panel A
Ax_handle[2] = Fig_handle.add_subplot(GS_handle[1:3, 8:10]) # Panel A
Ax_handle[3] = Fig_handle.add_subplot(GS_handle[1:3, 11:13]) # Panel A
Ax_handle[4] = Fig_handle.add_subplot(GS_handle[1:3, 14:16]) # Panel A
Ax_handle[5] = Fig_handle.add_subplot(GS_handle[4:5, 2:8]) # Panel B
Ax_handle[6] = Fig_handle.add_subplot(GS_handle[4:5, 10:16]) # Panel C

# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]
    
    Fig_handle.text(Text_X, Text_Y, Label, weight='bold',
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

add_panel_label(Ax_handle[0], 'A', (-0.1, 0.08))
add_panel_label(Ax_handle[5], 'B', (-0.1, 0.08))
add_panel_label(Ax_handle[6], 'C', (-0.08, 0.08))

#%% Collect data for the OD scatter plots

Panel_Scatter = {}
Panel_Scatter['All treatments'] = ['Arabinose 1E3', 'Arabinose 1E4', 'Arabinose 1E5', \
                        'Arabinose 1E6', 'Arabinose 1E7']
    
Panel_Scatter['All treatment labels'] = ['$10^3$ nM\nArabinose', '$10^4$ nM\nArabinose', '$10^5$ nM\nArabinose',
                        '$10^6$ nM\nArabinose', '$10^7$ nM\nArabinose']    
    
Panel_Scatter['Num Treatments'] = len(Panel_Scatter['All treatments'])

#%Calculate the raw correlations
Panel_Scatter['R values'] = np.zeros(Panel_Scatter['Num Treatments'])

for t, Treatment_Name in enumerate(Panel_Scatter['All treatments']):
    Panel_Scatter['R values'][t] = stats.pearsonr(Screen_Dict_384['Drug']['Wild-type control']['OD median'],
                                     Screen_Dict_384['Drug'][Treatment_Name]['OD median'])[0]

Panel_Scatter['R2 values'] = np.square(Panel_Scatter['R values'])


#%% Plot the OD scatter plots

Panel_Scatter['Axes'] = Ax_handle[0:5]

# Create a custom colormap for the scatter plots
Panel_Scatter['Base cmap'] = cm.get_cmap('cividis_r', 512)
Panel_Scatter['Custom cmap'] = ListedColormap(Panel_Scatter['Base cmap'](np.linspace(0.1, 1, 256)))

Panel_Scatter['Scatter kwargs'] = {'marker':'o', 'alpha':0.5, 's':1, 'cmap':Panel_Scatter['Custom cmap']}
Panel_Scatter['Text kwargs'] = {'horizontalalignment':'left', 'verticalalignment':'top',
                          'fontweight':'bold', 'color':KellyColors['Dark gray'], 'size':8, 'alpha':1}

# Make the plots
for Ax, Treatment_Name in zip(Panel_Scatter['Axes'], (Panel_Scatter['All treatments'])):
      
    # Choose the data sets to plot
    X_Data = Screen_Dict_384['Drug']['Wild-type control']['OD median']
    Y_Data = Screen_Dict_384['Drug'][Treatment_Name]['OD median']
    
    # Calculate the point density
    XY = np.vstack([X_Data,Y_Data])
    Z = stats.gaussian_kde(XY, bw_method=0.5)(XY)
    #Z = stats.gaussian_kde(XY)(XY)
    
    # Sort the points by density, so that the densest points are plotted on top
    idx = Z.argsort()
    X_Data, Y_Data, Z = X_Data[idx], Y_Data[idx], Z[idx]
    Ax.scatter(X_Data,Y_Data, c=Z, **Panel_Scatter['Scatter kwargs'])
    
# Label and style the axes
for Ax, Treatment_Name in zip(Panel_Scatter['Axes'], (Panel_Scatter['All treatments'])): 

    Ax.set_xlim([0, 0.45])
    Ax.set_xticks([0, 0.2, 0.4])
    Ax.set_ylim([0, 0.45])
    Ax.set_yticks([0, 0.2, 0.4])
    Ax.set_aspect('equal', 'box')
    Ax.grid(b=True, which='major', 
                      color=KellyColors['Medium gray'], alpha = 0.5)
    
    Ax.tick_params(axis='both', color='none')
    
# Annotate the plots
for t, Ax in enumerate(Panel_Scatter['Axes']):
    Ax.text(0.05, 0.95, '$r^2$ ' + str(round(Panel_Scatter['R2 values'][t], 2)),
           **Panel_Scatter['Text kwargs'], transform = Ax.transAxes,
           bbox=dict(facecolor='White', edgecolor='none',
                     boxstyle='Round, pad=0'))
    
    Ax.set_title(Panel_Scatter['All treatment labels'][t])

Ax_handle[2].text(0.5, -0.4, 'Wild-Type Control Growth (OD)', 
                  ha='center', transform=Ax_handle[2].transAxes)

Ax_handle[0].set_ylabel('TESEC ALR\nGrowth (OD)')

#%% Collect data for the SSMD violin plot

Panel_Violin = {};

Panel_Violin['All treatments'] = ['Wild-type control', 'Arabinose 1E3', 'Arabinose 1E4',
                        'Arabinose 1E5', 'Arabinose 1E6', 'Arabinose 1E7']

Panel_Violin['SSMD'] = [];  
for Treatment_Name in Panel_Violin['All treatments'] :
    Panel_Violin['SSMD'].append(Screen_Dict_384['Drug'][Treatment_Name]['SSMD'])
 
    
#%% Plot the SSMD violin plot

Ax = Ax_handle[5]

Panel_Violin['Color order'] = [KellyColors['Yellow green']] + 5 * [KellyColors['Purple']]

sns.violinplot(data=Panel_Violin['SSMD'], ax=Ax, 
                orient='v', linewidth=0, scale='area',
                bw=0.05, width=1,
                    palette=Panel_Violin['Color order'], saturation=1)

# Add marks for the DCS scores
Panel_Violin['DCS idx'] = 1130
Panel_Violin['DCS SSMD'] = [SSMD_Val[1130] for SSMD_Val in Panel_Violin['SSMD']]
Panel_Violin['DCS offset'] = np.arange(0,6) + 0.15
Ax.plot(Panel_Violin['DCS offset'],Panel_Violin['DCS SSMD'], 
        mec='black', marker='*', linestyle='none')


# Label and style the axes
Ax.set_xlabel('Arabinose (nM)')
Ax.set_ylabel('SSMD\nDMSO - Drug')

Ax.set_ylim([-18, 10])
Ax.set_yticks([-10, 0, 10])
Ax.set_xticklabels(('WTC', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'))

#%% Collect the sensitivity statistics

Panel_Sense = {}

Panel_Sense['All treatments'] = ['Wild-type control', 'Arabinose 1E3', 'Arabinose 1E4', 
                        'Arabinose 1E5', 'Arabinose 1E6', 'Arabinose 1E7', ]

Panel_Sense['True positive rate'] = np.array([
    Screen_Dict_384['Drug'][Treatment_Name]['True positive rate']
    for Treatment_Name in Panel_Sense['All treatments'] ])

Panel_Sense['False positive rate'] = np.array([
    Screen_Dict_384['Drug'][Treatment_Name]['False positive rate']
    for Treatment_Name in Panel_Sense['All treatments'] ])
    
Panel_Sense['True Negative Rate'] = np.array([
    Screen_Dict_384['Drug'][Treatment_Name]['True positive rate']
    for Treatment_Name in Panel_Sense['All treatments'] ])

Panel_Sense['False negative rate'] = np.array([
    Screen_Dict_384['Drug'][Treatment_Name]['False positive rate']
    for Treatment_Name in Panel_Sense['All treatments'] ])

Panel_Sense['False discovery rate'] = np.divide(
    Panel_Sense['False positive rate'], 
    (Panel_Sense['False positive rate'] + Panel_Sense['True positive rate']))

#%% Plot the sensitivity statistics

Ax = Ax_handle[6]

Ax.plot(range(1,6), Panel_Sense['True positive rate'][1:6], '-s', color=KellyColors['Orange'])
Ax.plot(range(1,6), Panel_Sense['False positive rate'][1:6], '-s', color=KellyColors['Dark green'])

Ax.plot(0, Panel_Sense['True positive rate'][0], 's', color=KellyColors['Orange'])
Ax.plot(0, Panel_Sense['False positive rate'][0],'s', color=KellyColors['Dark green'])

# Label and style the axes
Ax.set_yticks((0.25, 0.5, 0.75))
Ax.set_yticklabels(('25%', '50%', '75%'))
Ax.set_ylabel('Rate')

Ax.set_xticks(range(0,6))
Ax.set_xticklabels(('WTC', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'))
Ax.set_xlabel('Arabinose (nM)')

Ax.legend(labels=['True Positive', 'False Positive'],
    loc='upper left', frameon=False, ncol=3, 
          labelspacing=0, columnspacing=1, 
          handletextpad=0.5, borderaxespad=0,
          handlelength=1,
          bbox_to_anchor=(0.02, 0.9))

#%% Save the figure as a pdf

plt.savefig('Fig 3.pdf')