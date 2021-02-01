#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supporting information for:
    Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril
    
This script generates a figure in the manuscript from raw experimental data.

Figure 1 A TESEC strain for Mtb ALR shows differential senstivity to targeted inhibitors.
"""

#%% Import packages

# Data packages
from ast import literal_eval
import numpy as np
import pandas as pd
import json
from scipy.optimize import curve_fit

# Graphics packages
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator, FixedLocator

#%% Collect important metadata

Fig_meta = {};
Fig_meta['Data path'] = 'Source data/Fig 1/'
Fig_meta['Style path'] = 'Style files/'

# Load colors for plotting
with open(Fig_meta['Style path'] + 'Kelly colors.json') as json_file:
    KellyColors = json.load(json_file)

#%% Read data files used to make the plots

if not ('Variables_Loaded' in locals()):
    FCS_DF = pd.read_csv(Fig_meta['Data path'] + 'GFP-ALR single cell induction.csv', 
                         index_col=[0,1,2], converters={"FCS Data": literal_eval})

    # Read the OD data for the dose-response curve
    Dose_DF = pd.read_csv(Fig_meta['Data path'] + 'DCS dose response.csv', 
                          index_col=[0,1], header=[0,1])
    
    # Read the OD data for the DCS-arabinose curve
    IC50_DF = pd.read_csv(Fig_meta['Data path'] + 'DCS-arabinose IC50 curves.csv', 
                          index_col=[0, 1], header=[0,1])
    
    # Read the differential test
    DCS_Diff_DF = pd.read_csv(Fig_meta['Data path'] + 'DCS differential test.csv', 
                              index_col=[0], header=[0,1])
    
    Variables_Loaded = True
    
#%% Build the figure layout

# Apply the figure style sheet
plt.style.use(Fig_meta['Style path'] + 'Style Fig 1.mplstyle')

# Arrange the panels using GridSpec
Fig_handle = plt.figure(tight_layout=False)
Ax_handle = [None] * 6

Fig_meta['Rows'] = 10
Fig_meta['Cols'] = 18 
Fig_meta['Height ratios'] = [1] * Fig_meta['Rows'] 
Fig_meta['Width ratios'] = [1] * Fig_meta['Cols']
Fig_meta['Width ratios'][10] = 1
Fig_meta['Width ratios'][11] = 1

GS_handle = GridSpec(Fig_meta['Rows'], Fig_meta['Cols'], figure=Fig_handle, 
                     wspace=0, hspace = 0,
                     left=0, right=1, top=1, bottom=0,
                     height_ratios = Fig_meta['Height ratios'],
                     width_ratios = Fig_meta['Width ratios'])
                                     
Ax_handle[0] = Fig_handle.add_subplot(GS_handle[1:3, 2:4]) # Panel A
Ax_handle[1] = Fig_handle.add_subplot(GS_handle[1:3, 5:10]) # Panel B
Ax_handle[2] = Fig_handle.add_subplot(GS_handle[1:3, 12:17]) # Panel C
Ax_handle[3] = Fig_handle.add_subplot(GS_handle[6:8, 5:8]) # Panel D
Ax_handle[4] = Fig_handle.add_subplot(GS_handle[6:8, 10:12]) # Panel E
Ax_handle[5] = Fig_handle.add_subplot(GS_handle[5:8, 14:17]) # Panel F

# Label the panels using the upper-left points of their GridSpec positions
def add_panel_label(Panel, Label, Offsets):

    Text_X = Panel.get_position(Fig_handle).xmin + Offsets[0]
    Text_Y = Panel.get_position(Fig_handle).ymax + Offsets[1]
    
    Fig_handle.text(Text_X, Text_Y, Label, weight='bold', fontsize=8,
                    horizontalalignment='left', verticalalignment='top',
                    transform=Fig_handle.transFigure)

#add_panel_label(Ax_handle[0], 'A', (-0.06, 0.06))
add_panel_label(Ax_handle[1], 'B', (-0.06, 0.06))
add_panel_label(Ax_handle[2], 'C', (-0.06, 0.06))
add_panel_label(Ax_handle[3], 'D', (-0.06, 0.06))
add_panel_label(Ax_handle[4], 'E', (-0.06, 0.06))
add_panel_label(Ax_handle[5], 'F', (-0.06, 0.06))

#%% Panel A is a freehand drawing

Ax_handle[0].axis('off')

#%% Collect the data for the GFP violins

Panel_GFP = {}

# Collect the data used in the panel
Panel_GFP['Arabinose levels'] = np.log10(Dose_DF.columns.droplevel().values.astype(int))
Panel_GFP['FITC vals'] = FCS_DF.loc[('Rep A', 'FITC-A')]['FCS Data']

# Limit the plot data to the chosen range
Panel_GFP['X point min'] = 0
Panel_GFP['X point max'] = 23
Panel_GFP['X slice'] = slice(Panel_GFP['X point min'],Panel_GFP['X point max'])

Panel_GFP['Arabinose levels'] = Panel_GFP['Arabinose levels'][Panel_GFP['X slice']]
Panel_GFP['FITC vals'] = Panel_GFP['FITC vals'][Panel_GFP['X slice']]

Panel_GFP['Arabinose ticks'] = [2, 3, 4, 5, 6, 7, 8]
Panel_GFP['Arabinose tick labels'] = ['$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$', '$10^8$']

Panel_GFP['GFP limits'] = [3000,8500]
Panel_GFP['GFP ticks'] = [3000,8000]
Panel_GFP['GFP tick labels'] = ['3', '8']
Panel_GFP['GFP minor ticks']  = 1000

#%% Plot the GFP violins

Ax = Ax_handle[1]

Panel_GFP['Violin kwargs'] = {'showextrema':False, 'showmedians':True,
                              'bw_method':0.1, 'widths':0.15}

Panel_GFP['Violin handle'] = Ax.violinplot(Panel_GFP['FITC vals'], positions=Panel_GFP['Arabinose levels'],
                                          **Panel_GFP['Violin kwargs'])
Ax.set_ylabel('GFP (FITC)')
Ax.set_xticklabels([])

Ax.set_ylim(Panel_GFP['GFP limits'] )

# Recolor the violins and markers
for pc in Panel_GFP['Violin handle']['bodies']:
    pc.set_facecolor(KellyColors['Medium green'])
    pc.set_edgecolor('none')
    pc.set_alpha(1)

Panel_GFP['Violin handle']['cmedians'].set_color(KellyColors['Dark gray'])

Ax.set_xlabel('Arabinose (nM)')
Ax.set_xticks(Panel_GFP['Arabinose ticks'])
Ax.set_xticklabels(Panel_GFP['Arabinose tick labels'])
Ax.xaxis.set_minor_locator(MultipleLocator(0.25))

Ax.set_ylim(Panel_GFP['GFP limits'])
Ax.set_yticks(Panel_GFP['GFP ticks'])
Ax.set_yticklabels(Panel_GFP['GFP tick labels'] )
Ax.set_ylabel('GFP (FITC)')
Ax.yaxis.set_minor_locator(MultipleLocator(Panel_GFP['GFP minor ticks']))

#%% Collect the data for the differential signal curves

Panel_Signal = {}

# Copy parameters from GFP panel so they align
Panel_Signal['X slice'] = Panel_GFP['X slice']
Panel_Signal['Arabinose levels'] = Panel_GFP['Arabinose levels']

# Collect the data from the imported dataframes
Panel_Signal['With DCS mean'] = Dose_DF.loc['100 uM', 'Mean'].values
Panel_Signal['With DCS std'] = Dose_DF.loc['100 uM', 'Std'].values
Panel_Signal['No DCS mean'] = Dose_DF.loc['None', 'Mean'].values
Panel_Signal['No DCS std'] = Dose_DF.loc['None', 'Std'].values

# slice the parameters to their chosen range
Panel_Signal['Arabinose levels'] = Panel_Signal['Arabinose levels'][Panel_Signal['X slice']]
Panel_Signal['With DCS mean'] = Panel_Signal['With DCS mean'][Panel_Signal['X slice']]
Panel_Signal['With DCS std'] = Panel_Signal['With DCS std'][Panel_Signal['X slice']]
Panel_Signal['No DCS mean'] = Panel_Signal['No DCS mean'][Panel_Signal['X slice']]
Panel_Signal['No DCS std'] = Panel_Signal['No DCS std'][Panel_Signal['X slice']]

Panel_Signal['Mean difference'] = Panel_Signal['With DCS mean'] - Panel_Signal['No DCS mean'] 

Panel_Signal['Arabinose ticks'] = [2, 3, 4, 5, 6, 7, 8]
Panel_Signal['Arabinose tick labels'] = ['$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$', '$10^8$']

Panel_Signal['Growth limits'] = [0, 0.75]
Panel_Signal['Growth ticks'] = [0.2, 0.4]
Panel_Signal['Growth minor ticks'] = [0.1, 0.3]

#%% Plot the differential signal curves

Ax = Ax_handle[2]

Panel_Signal['Growth line kwargs'] = {'capsize':1, 'capthick':1}
Panel_Signal['Diff line kwargs'] = {'alpha':1}
Panel_Signal['Diff fill kwargs'] = {'interpolate':False, 'alpha':0.2}

# With DCS
Panel_Signal['Line with DCS'] = Ax.errorbar(Panel_Signal['Arabinose levels'], Panel_Signal['With DCS mean'] , yerr=Panel_Signal['With DCS std'],
            color=KellyColors['Medium blue'], **Panel_Signal['Growth line kwargs'])

# No DCS
Panel_Signal['Line no DCS'] = Ax.errorbar(Panel_Signal['Arabinose levels'], Panel_Signal['No DCS mean'] , yerr=Panel_Signal['No DCS std'], 
            color=KellyColors['Buff'], **Panel_Signal['Growth line kwargs'])

# Fill in the difference curve
Ax.fill_between(Panel_Signal['Arabinose levels'], Panel_Signal['With DCS mean'], 
                y2=Panel_Signal['No DCS mean'],
                where=(Panel_Signal['Mean difference'] < 0.08), 
                color=KellyColors['Buff'], **Panel_Signal['Diff fill kwargs'])

Ax.fill_between(Panel_Signal['Arabinose levels'], Panel_Signal['With DCS mean'], 
                y2=Panel_Signal['No DCS mean'],
                where=(Panel_Signal['Mean difference'] > 0), 
                color=KellyColors['Medium blue'], **Panel_Signal['Diff fill kwargs'])

# Label and style the axes
Ax.set_xlabel('Arabinose (nM)')
Ax.set_xticks(Panel_Signal['Arabinose ticks'])
Ax.set_xticklabels(Panel_Signal['Arabinose tick labels'])
Ax.xaxis.set_minor_locator(MultipleLocator(0.25))
Ax.set_ylabel('Growth (OD)')
Ax.set_ylim(Panel_Signal['Growth limits'] )
Ax.set_yticks(Panel_Signal['Growth ticks'])
Ax.yaxis.set_minor_locator(FixedLocator(Panel_Signal['Growth minor ticks']))

# Add the legend with the drug doses
Panel_Signal['Legend lines'] = (Panel_Signal['Line no DCS'], Panel_Signal['Line with DCS'])
Panel_Signal['Legend labels'] = ('No Drug', '0.1 mM DCS')
Ax.legend(Panel_Signal['Legend lines'], Panel_Signal['Legend labels'],
          loc='upper center', frameon=False, ncol=2, 
          labelspacing=0, columnspacing=1, 
          handletextpad=0.3, borderaxespad=0,
          handlelength=1, fontsize=7, markerscale=0.2,
          bbox_to_anchor=(0.5, 0.95))

# Annotate the high and low induction levels
Ax.annotate('', xy=(2, 0.05), xytext=(2,0.3),
            arrowprops={'arrowstyle':'->'})
Ax.annotate('Low\nArab.', xy=(2.3, 0.18), ha='left', va='center',
             bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.5)) 

Ax.annotate('', xy=(7, 0.05), xytext=(7,0.3),
            arrowprops={'arrowstyle':'->'})
Ax.annotate('High\nArab.', xy=(7.3, 0.18), ha='left', va='center',
            bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.5)) 

#%% Collect the data for the DCS-arabinose dose-response curves

Panel_DR = {}

# Stack up the dataframe into a 3d numpy array
Panel_DR['Dose Array'] = np.stack([IC50_DF.loc['Rep 1'].values, 
                          IC50_DF.loc['Rep 2'].values,
                          IC50_DF.loc['Rep 3'].values],
                          axis=-1)
                          
Panel_DR['Mean'] = np.mean(Panel_DR['Dose Array'], axis=2)
Panel_DR['Std'] = np.std(Panel_DR['Dose Array'], axis=2)
Panel_DR['DCS (mM)'] = IC50_DF.columns.droplevel().values.astype(float)
Panel_DR['Arabinose (nM)'] = IC50_DF.loc['Rep 1'].index.values.astype(float)
Panel_DR['Log DCS'] = np.log10(Panel_DR['DCS (mM)'])
Panel_DR['Log DCS Smooth'] = np.arange(min(Panel_DR['Log DCS']),
                                                 max(Panel_DR['Log DCS']),
                                                 0.05)
Panel_DR['Log arabinose'] = np.log10(Panel_DR['Arabinose (nM)'])

Panel_DR['DCS limits'] = np.log10([5, 2000])
Panel_DR['DCS ticks'] = np.log10([20, 200, 2000])
Panel_DR['DCS tick labels'] = ['0.02', '0.2', '2']

Panel_DR['Arabinose ticks'] = [3, 5, 7]
Panel_DR['Arabinose tick labels'] = ['$10^3$','$10^5$','$10^7$']

Panel_DR['Arabinose Legend Labels'] = ['$10^2$','$10^3$','$10^4$','$10^5$','$10^6$','$10^7$','$10^8$']

Panel_DR['Growth limits'] = [0, 0.8]
Panel_DR['Growth ticks'] = [0, 0.3, 0.6]
Panel_DR['Growth minor ticks'] = 0.1

#%% Plot the DCS-arabinose dose-response curves

Ax = Ax_handle[3]

# Define the colors and properties of the lines
Panel_DR['Colormap'] = plt.cm.inferno_r
Panel_DR['Num lines'] = 8
Panel_DR['Line colors'] = Panel_DR['Colormap'](np.linspace(0, 1, Panel_DR['Num lines']))
Ax.set_prop_cycle(color=Panel_DR['Line colors'])

# Create the lines
for r in range(8):
    Ax.errorbar(Panel_DR['Log DCS'], Panel_DR['Mean'][r, :], yerr=Panel_DR['Std'][r, :])
    
# Label and style the axes
Ax.set_xticks(Panel_DR['DCS ticks'])
Ax.set_xticklabels(Panel_DR['DCS tick labels'])
Ax.set_xlim(Panel_DR['DCS limits'])
Ax.set_xlabel('DCS (mM)')

Ax.set_yticks(Panel_DR['Growth ticks'])
Ax.set_ylim(Panel_DR['Growth limits'])
Ax.set_ylabel('Growth (OD)')

# Create an axis for the colorbar
Panel_DR['Color norm']  = mpl.colors.Normalize(vmin=2, vmax=8)
Panel_DR['Colorbar ax'] = inset_axes(Ax,
                    width="100%",
                    height="100%",
                    bbox_to_anchor=(0.1, 1.2, 0.8, 0.1),
                    bbox_transform=Ax.transAxes,
                    borderpad=0,
                    )
Panel_DR['Colorbar handle'] = mpl.colorbar.ColorbarBase(
    Panel_DR['Colorbar ax'],
    cmap=Panel_DR['Colormap'],
    norm=Panel_DR['Color norm'],
    orientation='horizontal')
Panel_DR['Colorbar handle'].set_ticks([2,3,4,5,6,7,8])
Panel_DR['Colorbar handle'].set_ticklabels([])
Ax.yaxis.set_minor_locator(MultipleLocator(Panel_DR['Growth minor ticks']))
Panel_DR['Colorbar handle'].ax.xaxis.set_label_position('bottom')
Panel_DR['Colorbar handle'].ax.set_xlabel('Arabinose (nM)')

# Mark the high and low induction levels
Ax.annotate('$10^2$', xy=(0, 1.25), xycoords = 'axes fraction',
            ha='center', va='center')

Ax.annotate('$10^8$', xy=(1, 1.25), xycoords = 'axes fraction',
            ha='center', va='center')

#%% Fit the dose-response curves to a logistic function

Panel_IC50 = {}

def logistic_func(x, Lmin, Lmax, k, IC50):
    return Lmin + (Lmax / (1 + np.exp(-k *(x-IC50))))

Cfit = {}
Cfit['Num fits'] = 8
Cfit['Bounds'] = ([0, 0, -10, -5], [0.01, 1, 10, 5])
Cfit['Fit params'] = 8 * [object()]
Cfit['Param cov'] = 8 * [object()]

for i in range(Cfit['Num fits']):
    Cfit['Fit params'][i], Cfit['Param cov'][i] = curve_fit(logistic_func, 
        Panel_DR['Log DCS'] , Panel_DR['Mean'][i, :], bounds=Cfit['Bounds'])

Panel_IC50['Fit IC50'] = [None] * Cfit['Num fits']
Panel_IC50['Fit IC50 95 CI'] = [None] * Cfit['Num fits']

# Collect the IC50 fits and confidence intervals for each arabinose level
for i in range(Cfit['Num fits']):
    Panel_IC50['Fit IC50'][i] = Cfit['Fit params'][i][3]
    Panel_IC50['Fit IC50 95 CI'][i] = 1.96 * np.sqrt(np.diag(Cfit['Param cov'][i])[3])
    
#%% Collect the data for the IC50 curves

Panel_IC50['Arabinose (nM)'] = Panel_DR['Arabinose (nM)']
Panel_IC50['Log arabinose'] = np.log10(Panel_IC50['Arabinose (nM)'])

Panel_IC50['Arabinose limits'] = [0.5, 8.5]
Panel_IC50['Arabinose ticks'] = [2, 5, 8]
Panel_IC50['Arabinose tick labels'] = ['$10^2$', '$10^5$', '$10^8$']

Panel_IC50['DCS limits'] = np.log10([5, 2000])
Panel_IC50['DCS ticks'] = np.log10([20, 200, 2000])
Panel_IC50['DCS tick labels'] = ['0.02', '0.2', '2']

Panel_IC50['Error bar kwargs'] = {'color':KellyColors['Violet']}

#%% Plot the IC50 curves

Ax = Ax_handle[4]

Ax.errorbar(Panel_IC50['Log arabinose'], Panel_IC50['Fit IC50'], yerr=Panel_IC50['Fit IC50 95 CI'], 
            **Panel_IC50['Error bar kwargs'])

# Label and style the axes
Ax.set_xlim(Panel_IC50['Arabinose limits'])
Ax.set_xticks(Panel_IC50['Arabinose ticks'])
Ax.set_xticklabels(Panel_IC50['Arabinose tick labels'])
Ax.xaxis.set_minor_locator(MultipleLocator(1))
Ax.set_xlabel('Arabinose (nM)')
Ax.set_ylim(Panel_IC50['DCS limits'])
Ax.set_yticks(Panel_IC50['DCS ticks'])
Ax.set_yticklabels(Panel_IC50['DCS tick labels'])
Ax.set_ylabel('DCS IC$_{50}$ (mM)')

#%% Collect the data for the differential screen proof-of-concept

Panel_POC = {}
Panel_POC['High ODs'] = IC50_DF.xs(1E7, axis=0, level=1, drop_level=False)
Panel_POC['High mean'] = np.mean(Panel_POC['High ODs']).values
Panel_POC['High 95 conf'] = (1.96/np.sqrt(3)) * np.std(Panel_POC['High ODs'].values, axis=0) 

Panel_POC['Low ODs'] = IC50_DF.xs(1E2, axis=0, level=1, drop_level=False)
Panel_POC['Low mean'] = np.mean(Panel_POC['Low ODs']).values
Panel_POC['Low 95 conf'] = (1.96/np.sqrt(3)) * np.std(Panel_POC['Low ODs'].values, axis=0) 

Panel_POC['DCS (μM)'] = np.log10(Panel_POC['High ODs'].columns.get_level_values(1).values.astype(float))

Panel_POC['DCS ticks'] = np.log10([20, 200, 2000])
Panel_POC['DCS tick labels'] = ['0.02', '0.2', '2']

Panel_POC['Growth limits'] = [-0.05,0.75]
Panel_POC['Growth ticks'] = [0, 0.3, 0.6]
Panel_POC['Growth minor ticks'] = 0.1

#%% Plot the differential screen proof-of-concept

Ax = Ax_handle[5]

Panel_POC['Diff kwargs'] = {'marker':'o', 'cmap':'viridis_r', 
                             'linewidth':0.5, 'edgecolors':KellyColors['Dark gray']}

Panel_POC['Scatter handle'] = Ax.scatter(Panel_POC['Low mean'], Panel_POC['High mean'], 
                                          c=Panel_POC['DCS (μM)'], **Panel_POC['Diff kwargs'])

# Create an axis for the colorbar
Panel_POC['Colorbar ax'] = inset_axes(Ax,
                    width="100%",
                    height="100%",
                    bbox_to_anchor=(1.16, 0.1, 0.08, 0.8),
                    bbox_transform=Ax.transAxes,
                    borderpad=0,
                    )

Panel_POC['Colorbar handle'] = plt.colorbar(Panel_POC['Scatter handle'], 
                                              cax=Panel_POC['Colorbar ax'],
                                              orientation='vertical',
                                              anchor=(1, 0.5))


Panel_POC['Colorbar handle'].set_ticks(Panel_POC['DCS ticks'])
Panel_POC['Colorbar handle'].set_ticklabels([])
Panel_POC['Colorbar handle'].ax.yaxis.set_label_position('left')
Panel_POC['Colorbar handle'].ax.yaxis.set_ticks_position('right')
Panel_POC['Colorbar handle'].ax.set_ylabel('DCS (mM)')
Panel_POC['Colorbar handle'].ax.invert_yaxis()

# Label and style the axes
Ax.set_xlabel('Low Arab. (OD)')
Ax.set_ylabel('High Arab. (OD)')
Ax.set_xlim(Panel_POC['Growth limits'])
Ax.set_xticks(Panel_POC['Growth ticks'])
Ax.xaxis.set_minor_locator(MultipleLocator(Panel_POC['Growth minor ticks']))
Ax.set_ylim(Panel_POC['Growth limits'])
Ax.set_yticks(Panel_POC['Growth ticks'])
Ax.yaxis.set_minor_locator(MultipleLocator(Panel_POC['Growth minor ticks']))
Ax.set_aspect('equal')

Ax.annotate('0.02', xy=(1.2, 1), xycoords = 'axes fraction',
            ha='center', va='center')

Ax.annotate('2', xy=(1.2, 0), xycoords = 'axes fraction',
            ha='center', va='center')

#%% Save the figure as a pdf

plt.savefig('Fig 1.pdf')
