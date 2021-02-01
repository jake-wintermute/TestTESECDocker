#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril
    
This script generates a figure in the manuscript from raw experimental data.

Figure 2 A screen for targeted Mtb-ALR inhibitors identifies benazepril. 
"""

#%% Import packages

# Data packages
import numpy as np
import pickle
import json

# Graphics packages
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 2/'
Fig_meta['Style path'] = 'Style files/'

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)

#%% Read data files used to make the plots

with open(Fig_meta['Data path'] + 'TESEC ALR Screen Dict 96.pickle', 'rb') as Pickle_File_Handle:
    Screen_Dict_96  = pickle.load(Pickle_File_Handle)
    
with open(Fig_meta['Data path'] + 'TESEC ALR Screen Dict 384.pickle', 'rb') as Pickle_File_Handle:
    Screen_Dict_384  = pickle.load(Pickle_File_Handle)

#%% Build the figure layout

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] + 'Style Fig 2.mplstyle')

# Arrange the panels using GridSpec
Fig_handle = plt.figure(tight_layout=False)
Ax_handle = [None] * 5

Fig_meta['Rows'] = 7
Fig_meta['Cols'] = 11 
Fig_meta['Height ratios'] = [1,1,1,1,1.5,1.5,1]
Fig_meta['Width ratios'] = [1,1,1,1,1,1,1,1,1,1,0.5]

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle, 
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'])
                                     
Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1:3, 1:3]) # Panel A
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[1:3, 4:6]) # Panel B
Ax_handle[2] = Fig_handle.add_subplot(GS_handle[2, 7:10]) # Panel C
Ax_handle[3] = Fig_handle.add_subplot(GS_handle[5, 1:5]) # Panel D
Ax_handle[4] = Fig_handle.add_subplot(GS_handle[5, 6:10]) # Panel E
Legend_position = GS_handle[5, 10]

# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]
    
    Fig_handle.text(Text_X, Text_Y, Label, weight='bold', fontsize=8,
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

add_panel_label(Ax_handle[0], 'A', (-0.08, 0.08))
add_panel_label(Ax_handle[1], 'B', (-0.08, 0.08))
add_panel_label(GS_handle[1, 7], 'C', (-0.08, 0.08))
add_panel_label(Ax_handle[3], 'D', (-0.08, 0.08))
add_panel_label(Ax_handle[4], 'E', (-0.08, 0.08))

#%% Select the hit compounds

treatment_1 = 'Low arabinose'
treatment_2 = 'High arabinose'

Thresholds = {}
Thresholds['OD low'] = 0.1
Thresholds['OD high'] = 0.15
Thresholds['SSMD'] = 5

Hit_Selector = {}
Hit_Selector['Boolean'] = {}

Hit_Selector['Boolean']['Killed with treatment 1'] = \
    Screen_Dict_96['Drug'][treatment_1]['OD mean'] < Thresholds['OD low']
    
Hit_Selector['Boolean']['Survives with treatment 2'] = \
    Screen_Dict_96['Drug'][treatment_2]['OD mean'] > Thresholds['OD high']   

Hit_Selector['Boolean']['Significant SSMD'] = \
    Screen_Dict_96['Pairwise'][treatment_2][treatment_1]['SSMD'] > Thresholds['SSMD']

Hit_Selector['Boolean']['Hits'] = \
    Hit_Selector['Boolean']['Significant SSMD'] & \
    Hit_Selector['Boolean']['Killed with treatment 1'] & \
    Hit_Selector['Boolean']['Survives with treatment 2']

# Collect the library info for the hit compounds
Hit_Selector['Drug index'] = {}
Hit_Selector['Drug index']['Hits'] = np.where(Hit_Selector['Boolean']['Hits'])[0]
Hit_Selector['Drug index']['Non-hits'] = list(set(range(len(Screen_Dict_96['Library info']))) - \
                                        set(Hit_Selector['Drug index']['Hits']))

Hit_Selector['Profile'] = Screen_Dict_96['Library info'].iloc[Hit_Selector['Drug index']['Hits']]
Hit_Selector['Names'] = Screen_Dict_96['Library info']['Chemical name'].iloc[Hit_Selector['Drug index']['Hits']].to_list()

#%% Build a dict of all the hit compound properties

# Collect the screen data for all the hits
Hit_Dict = {}
Hit_Dict['All drug names'] = Hit_Selector['Names']
for Drug_Name in Hit_Dict['All drug names']:
    Hit_Dict[Drug_Name] = {}
    Hit_Dict[Drug_Name]['Drug index'] = np.where(Screen_Dict_96['Library info']['Chemical name'] == Drug_Name)[0][0]
    Hit_Dict[Drug_Name]['OD low mean'] = Screen_Dict_96['Drug'][treatment_1]['OD mean'][Hit_Dict[Drug_Name]['Drug index']]
    Hit_Dict[Drug_Name]['OD low std'] = Screen_Dict_96['Drug'][treatment_1]['OD std'][Hit_Dict[Drug_Name]['Drug index']]
    Hit_Dict[Drug_Name]['OD high mean'] = Screen_Dict_96['Drug'][treatment_2]['OD mean'][Hit_Dict[Drug_Name]['Drug index']]
    Hit_Dict[Drug_Name]['OD high std'] = Screen_Dict_96['Drug'][treatment_2]['OD std'][Hit_Dict[Drug_Name]['Drug index']]
    Hit_Dict[Drug_Name]['OD diff'] = Screen_Dict_96['Pairwise'][treatment_2][treatment_1]['OD mean diff'][Hit_Dict[Drug_Name]['Drug index']]
    Hit_Dict[Drug_Name]['SSMD'] = Screen_Dict_96['Pairwise'][treatment_2][treatment_1]['SSMD'][Hit_Dict[Drug_Name]['Drug index']]

# Set the base graphical features for plotting all hits. Interesting hits will get updated features
for Drug_Name in Hit_Dict['All drug names']:
    Hit_Dict[Drug_Name]['Short name'] = Drug_Name.split()[0]
    Hit_Dict[Drug_Name]['Face color'] = KellyColors['Buff']
    Hit_Dict[Drug_Name]['Edge color'] = 'none'
    Hit_Dict[Drug_Name]['Edge width'] = 0
    Hit_Dict[Drug_Name]['Bar color'] = KellyColors['Buff']
    Hit_Dict[Drug_Name]['Line color'] = KellyColors['Buff']
    Hit_Dict[Drug_Name]['Marker'] = 'o'
    Hit_Dict[Drug_Name]['Size'] = 6

# Set the special properties for the interesting hits
Drug_Name = 'D-cycloserine'
Hit_Dict[Drug_Name]['Short name'] = 'DCS'
Hit_Dict[Drug_Name]['Face color'] = 'none'
Hit_Dict[Drug_Name]['Edge color'] = KellyColors['Medium blue']
Hit_Dict[Drug_Name]['Edge width'] = 2
Hit_Dict[Drug_Name]['Bar color'] = KellyColors['Medium blue']
Hit_Dict[Drug_Name]['Line color'] = KellyColors['Medium blue']
Hit_Dict[Drug_Name]['Marker'] = '+'
Hit_Dict[Drug_Name]['Size'] = 6
    
Drug_Name = 'Amlexanox'
Hit_Dict[Drug_Name]['Face color'] = 'none'
Hit_Dict[Drug_Name]['Edge color'] = KellyColors['Orange yellow']
Hit_Dict[Drug_Name]['Edge width'] = 2
Hit_Dict[Drug_Name]['Bar color'] = KellyColors['Orange yellow']
Hit_Dict[Drug_Name]['Line color'] = KellyColors['Orange yellow']
Hit_Dict[Drug_Name]['Marker'] = 'o'
Hit_Dict[Drug_Name]['Size'] = 6

Drug_Name = 'Benazepril HCl'
Hit_Dict[Drug_Name]['Face color'] = 'none'
Hit_Dict[Drug_Name]['Edge color'] = KellyColors['Dark pink']
Hit_Dict[Drug_Name]['Edge width'] = 2
Hit_Dict[Drug_Name]['Bar color'] = KellyColors['Dark pink']
Hit_Dict[Drug_Name]['Line color'] = KellyColors['Dark pink']
Hit_Dict[Drug_Name]['Marker'] = 'o'
Hit_Dict[Drug_Name]['Size'] = 6

# Group the drugs
Hit_Dict['Positive control'] = ('D-cycloserine',)
Hit_Dict['Unknown mechanism'] = ('Amlexanox', 'Benazepril HCl')
Hit_Dict['Beta lactams'] = 	('Cefotiam hydrochloride', 'Aztreonam', 'Ceforanide', 'Loracarbef',
                             'Cefuroxime axetil', 'Cefpodoxime proxetil', 'Ceftibuten')

#%% Collect data for the OD scatter plot showing the screen hits

Panel_Screen = {}
Panel_Screen['OD 1 all'] = Screen_Dict_96['Drug'][treatment_1]['OD600 Median']
Panel_Screen['OD 2 all'] = Screen_Dict_96['Drug'][treatment_2]['OD600 Median']

#Zero out the negative ODs for plotting
Negative_OD_idx = Panel_Screen['OD 1 all'] < 0
Panel_Screen['OD 1 all'][Negative_OD_idx] = 0
Negative_OD_idx = Panel_Screen['OD 2 all'] < 0
Panel_Screen['OD 2 all'][Negative_OD_idx] = 0

Panel_Screen['OD 1 non-hits'] = Panel_Screen['OD 1 all'][Hit_Selector['Drug index']['Non-hits']]
Panel_Screen['OD 2 non-hits'] = Panel_Screen['OD 2 all'][Hit_Selector['Drug index']['Non-hits']]

#%% Plot the OD scatter plot showing the screen hits

Ax = Ax_handle[0]

# Plot the non-hits
Ax.scatter(Panel_Screen['OD 1 non-hits'], Panel_Screen['OD 2 non-hits'],
                    c=KellyColors['Medium gray'], edgecolors='none',
                          alpha=0.3)

# Plot the hits 
for Drug_Name in (Hit_Dict['Beta lactams'] + 
                  Hit_Dict['Unknown mechanism'] +
                  Hit_Dict['Positive control']):
    
    Marker_Properties = {'mfc' : Hit_Dict[Drug_Name]['Face color'],
                         'mec' : Hit_Dict[Drug_Name]['Edge color'], 
                         'marker' : Hit_Dict[Drug_Name]['Marker'],
                         'ms' : Hit_Dict[Drug_Name]['Size'],
                         'markeredgewidth' : Hit_Dict[Drug_Name]['Edge width'], 
                         'alpha' : 1}
    
    Ax.plot(Hit_Dict[Drug_Name]['OD low mean'],
            Hit_Dict[Drug_Name]['OD high mean'],
               **Marker_Properties)

# Label and style the axes
Ax.set_xlabel('Low Arab. (OD)')
Ax.set_xticks((0, 0.2,0.4))
Ax.set_ylabel('High Arab. (OD)')
Ax.set_yticks((0, 0.2,0.4))
Ax.set_aspect('equal')

#%% Collect data for the SSMD scatter plot

Panel_SSMD = {}
Panel_SSMD['OD diff all'] = Screen_Dict_96['Pairwise'][treatment_2][treatment_1]['OD mean diff']
Panel_SSMD['OD diff hits'] = Panel_SSMD['OD diff all'][Hit_Selector['Drug index']['Hits']]
Panel_SSMD['OD diff non-hits'] = Panel_SSMD['OD diff all'][Hit_Selector['Drug index']['Non-hits']]
Panel_SSMD['SSMD all'] = abs(Screen_Dict_96['Pairwise'][treatment_2][treatment_1]['SSMD'])
Panel_SSMD['SSMD hits'] = Panel_SSMD['SSMD all'][Hit_Selector['Drug index']['Hits']]
Panel_SSMD['SSMD non-hits'] = Panel_SSMD['SSMD all'][Hit_Selector['Drug index']['Non-hits']]

#%% Plot the SSMD scatter plot

Ax = Ax_handle[1]

# Plot the non-hits
Ax.scatter(Panel_SSMD['OD diff non-hits'], Panel_SSMD['SSMD non-hits'],
                          c=KellyColors['Medium gray'], edgecolors='none',
                          alpha=0.3)
    
# Plot the hits 
for Drug_Name in (Hit_Dict['Beta lactams'] + 
                  Hit_Dict['Positive control'] +
                  Hit_Dict['Unknown mechanism']):

    Marker_Properties = {'mfc' : Hit_Dict[Drug_Name]['Face color'],
                         'mec' : Hit_Dict[Drug_Name]['Edge color'], 
                         'marker' : Hit_Dict[Drug_Name]['Marker'],
                         'ms' : Hit_Dict[Drug_Name]['Size'],
                         'markeredgewidth' : Hit_Dict[Drug_Name]['Edge width']} 
    
    Ax.plot(Hit_Dict[Drug_Name]['OD diff'], 
            Hit_Dict[Drug_Name]['SSMD'],
               **Marker_Properties)

#Label and style the axes
Ax.set_xlabel('Growth Difference\n(OD$_{HIGH}$ - OD$_{LOW}$)')
Ax.set_ylabel('|SSMD|\nOD$_{HIGH}$ - OD$_{LOW}$')
Ax.set_ylim((0, 50))


#%% Collect data for the bar plot of hit compounds

Panel_Bars = {}
Panel_Bars['OD 1 means'] = Screen_Dict_96['Drug unsubtracted'][treatment_1]['OD mean'][Hit_Selector['Drug index']['Hits']]
Panel_Bars['OD 1 stds'] = Screen_Dict_96['Drug unsubtracted'][treatment_1]['OD std'][Hit_Selector['Drug index']['Hits']]
Panel_Bars['OD 2 Means'] = Screen_Dict_96['Drug unsubtracted'][treatment_2]['OD mean'][Hit_Selector['Drug index']['Hits']]
Panel_Bars['OD 2 Stds'] = Screen_Dict_96['Drug unsubtracted'][treatment_2]['OD std'][Hit_Selector['Drug index']['Hits']]

#%% Plot the bar plot of hit compounds

Ax = Ax_handle[2]

Panel_Bars = {}

Panel_Bars['Bar thickness'] = 0.4;
Panel_Bars['Bar offset'] = 0.4;
Panel_Bars['Bar positions'] = (1,) + (3, 4, 5, 6, 7, 8, 9) + (11, 12)
Panel_Bars['Hit drugs']  =  Hit_Dict['Positive control'] + Hit_Dict['Beta lactams'] + Hit_Dict['Unknown mechanism']

Panel_Bars['Drug labels'] = [Hit_Dict[Drug_Name]['Short name'] for Drug_Name in Panel_Bars['Hit drugs']]

for Bar_Pos, Drug_Name in zip(Panel_Bars['Bar positions'], Panel_Bars['Hit drugs']):
    
    # Plot the high induction ODs
    Ax.bar(Bar_Pos, Hit_Dict[Drug_Name]['OD high mean'],
              yerr = Hit_Dict[Drug_Name]['OD high std'],
              width = Panel_Bars['Bar thickness'], label = Hit_Dict[Drug_Name]['Short name'], 
              color = Hit_Dict[Drug_Name]['Bar color'], capsize=1)
    
    # Plot the low induction ODs
    Ax.bar(Bar_Pos+Panel_Bars['Bar offset'], Hit_Dict[Drug_Name]['OD low mean'],
          yerr = Hit_Dict[Drug_Name]['OD low std'],
          width = Panel_Bars['Bar thickness'],
          color = KellyColors['Medium gray'], capsize=1)
    
# Label and style the axes
Ax.set_xticks(np.array(Panel_Bars['Bar positions'] ) + (Panel_Bars['Bar offset'] /2))
Ax.set_xticklabels(labels = Panel_Bars['Drug labels'], rotation=45, ha='right')
Ax.set_xlim([0.4, 12.8])
Ax.set_ylabel('Growth (OD)')

#%% Collect data for the chem-gen growth traces

Panel_Trace = {};

Panel_Trace['All trace treatments'] = ['Arabinose 1E3', 'Arabinose 1E4', 'Arabinose 1E5',
                        'Arabinose 1E6', 'Arabinose 1E7']
     
Panel_Trace['OD mean'] = np.vstack([Screen_Dict_384['Drug'][Trace_Treatment]['OD mean']
                          for Trace_Treatment in Panel_Trace['All trace treatments']])
    
Panel_Trace['OD std'] = np.vstack([Screen_Dict_384['Drug'][Trace_Treatment]['OD std']
                          for Trace_Treatment in Panel_Trace['All trace treatments']])

#%% Plot the chem-gen growth traces

Ax = Ax_handle[3]

# Plot all the traces
Ax.plot(Panel_Trace['OD mean'],
                  color=KellyColors['Dark gray'], alpha=0.02)

# Label and style the axes
Ax.set_ylim((0,0.45))
Ax.set_yticks((0, 0.2, 0.4))
Ax.set_ylabel('Growth (OD)')
Ax.set_xticks([0, 1, 2, 3, 4])
Ax.set_xticklabels(('$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'))
Ax.set_xlabel('Arabinose (nM)')

#%% Collect data for the chem-gen trace highlights

Panel_Hits = {}

Panel_Hits['OD mean'] = np.vstack([Screen_Dict_384['Drug'][Trace_Treatment]['OD mean']
                          for Trace_Treatment in Panel_Trace['All trace treatments']])
    
Panel_Hits['OD std'] = np.vstack([Screen_Dict_384['Drug'][Trace_Treatment]['OD std']
                          for Trace_Treatment in Panel_Trace['All trace treatments']])

Panel_Hits['Highlights'] = {};
Panel_Hits['Highlights']['Drug names'] = ['Amlexanox', 'D-cycloserine', 'Benazepril HCl']

Panel_Hits['Highlights']['Indices'] = []
for Drug_Name in Panel_Hits['Highlights']['Drug names']:
    Drug_Idx = np.where(Screen_Dict_384['Library info']['Chemical name'].str.contains(Drug_Name, regex=False))[0][0]
    Panel_Hits['Highlights']['Indices'].append(Drug_Idx)

for (Drug_Idx, Drug_Name) in zip(Panel_Hits['Highlights']['Indices'], Panel_Hits['Highlights']['Drug names']):
    Panel_Hits['Highlights'][Drug_Name] = {};
    Panel_Hits['Highlights'][Drug_Name]['OD mean'] = Panel_Hits['OD mean'][:, Drug_Idx]
    Panel_Hits['Highlights'][Drug_Name]['OD std'] = Panel_Hits['OD std'][:, Drug_Idx]


#%% Plot the chem-gen trace highlights

Ax = Ax_handle[4]

# Plot the interesting traces
for Drug_Name in Panel_Hits['Highlights']['Drug names']:
    Ax.errorbar(range(5), Panel_Hits['Highlights'][Drug_Name]['OD mean'],
                      yerr = Panel_Hits['Highlights'][Drug_Name]['OD std'],
                      label = Hit_Dict[Drug_Name]['Short name'], 
                      color = Hit_Dict[Drug_Name]['Line color'], linewidth=2)

# Label and style the axes
Ax.set_xticks([0, 1, 2, 3, 4])
Ax.set_xticklabels(('$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'))
Ax.set_xlabel('Arabinose (nM)')
Ax.set_ylim((0,0.45))
Ax.set_yticks((0, 0.2, 0.4))
Ax.set_ylabel('Growth (OD)')

Ax.legend(loc='upper center', frameon=False, ncol=3, 
          labelspacing=0, columnspacing=1, 
          handletextpad=0.3, borderaxespad=0,
          handlelength=1,
          bbox_to_anchor=(0.5, 1))


#%% Save the figure as a pdf

plt.savefig('Fig 2.pdf')