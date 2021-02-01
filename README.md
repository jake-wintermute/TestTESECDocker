# Supporting data for: Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril

## Abbreviations used in these files
TESEC: Target Essential Surrogate E. coli \
Mtb: Mycobacterium tuberculosis \
ALR: Alanine racemase \
DCS: D-Cycloserine \
Dala: D-alanine \
IC50: Half-maximal inhibitory concentration \
GFP: Green Fluorescent Protein

## Background & summary of the project
The associated reference describes our work developing TESEC Mtb ALR, a genetically engineered strain of E. coli expressing the enzyme ALR derived from Mtb. We used the TESEC Mtb ALR strain in a high-throughput drug screen and identified benazepril as targeted inhibitor of the ALR enzyme. We then performed additional experiments to characterize the activity of benazepril against E. coli, Mtb and purified enzymes.

These files include growth measurements, biochemical assays and other forms of biological data. They are packaged together with scripts used to analyze the data and present them in figures. Our goal in creating this archive was to present our complete analysis pipeline in the spirit of open science. It is not intended to be a stand-alone resource - the associated reference provides protocols, units of measurement and other essential technical context.

Our scripts were written for Python 3.8. The raw data is presented as human-readable .csv files intended to be imported as Pandas DataFrames. Some data is also packaged as Python dictionaries saved with the Pickle package.

## References
Low-cost drug discovery with engineered E. coli reveals an anti-mycobacterial activity of benazepril

## Authors
Nadine Bongaerts¹, Zainab Edoo², Ayan A. Abukar¹, Xiaohu Song¹, Sebastian Sosa Carrillo³, Ariel B. Lindner¹ & Edwin H Wintermute¹

¹ INSERM U1284, Université de Paris, Center for Research and Interdisciplinarity (CRI), 75004 Paris, France \
² Centre de Recherche des Cordeliers, CRC, 75006, Paris, France \
³ USR 3756 IP CNRS, Institut Pasteur, Paris, France

## How to use this data
These python scripts can be run with Python 3.8.5 to generate the figures in our manuscript from raw data. They can be run in the python environment of your choice with access to standard packages, for example matplotlib, listed in the file requirements.txt.

Alternately, the scripts can be run using on any system with Python 3 installed using the supplied virtual environment, pyenv.pex. To run these scripts from the command line type for example:

python ./pyenv.pex GenerateFig1.py

The script will use the raw data in the 'Source Data' folder to generate and save a .pdf version of the figure.

## Files included in this dataset

### Figures generated from the script
Fig 1.pdf, Fig 2.pdf, Fig 3.pdf, Fig 4.pdf, Fig 5.pdf, Fig 6.pdf

### Raw data files used to create the figures
DCS differential test.csv \
DCS dose response.csv \
DCS-arabinose IC50 curves.csv \
GFP-ALR single cell induction.csv \
Prestwick Info 96 Format.csv \
Prestwick Info 384 Format.csv \
TESEC ALR Screen DF 96.csv \
TESEC ALR Screen DF 384.csv \
TESEC ALR Screen Dict 96.pickle \
TESEC ALR Screen Dict 384.pickle \
TESEC ChemGen Profiles.csv \
TESEC Dala Rescue.csv \
TESEC Efflux Effect.csv \
Msmeg Dala Rescue.csv \
Msmeg IC50 Curves.csv \
In Vitro Hanes Woolf.csv \
In Vitro Michaelis Menten.csv

### Python scripts used to create the figures
GenerateFig1.py, GenerateFig2.py, GenerateFig3.py, GenerateFig4.py, GenerateFig5.py, GenerateFig6.py

### Python environment details
pyenv.pex \
requirements.txt

### Figure styling tools
Style Fig 1.mplstyle, Style Fig 2.mplstyle, Style Fig 3.mplstyle, Style Fig 4.mplstyle, Style Fig 5.mplstyle, Style Fig 6.mplstyle \
Kellycolors.json

## Description of file contents

Files according to their use in generating figures of the associated manuscript.

### Figure 1 A TESEC strain for Mtb ALR shows differential sensitivity to targeted inhibitors

**DCS differential test.csv** \
Growth data for the TESEC Mtb ALR strain with a range of DCS levels and ALR induction at either high or low levels. Used to demonstrate that TESEC strains show differential sensitivity to a targeted inhibitor under screen-like conditions.

**DCS dose response.csv** \
Growth data for the TESEC Mtb ALR strain for a range of arabinose levels, with and without DCS treatment. Used to identify the concentrations of arabinose that produce the largest differential response.

**DCS-arabinose IC50 curves.csv** \
Growth data for the TESEC Mtb ALR strain varying both DCS and arabinose. Used to fit the IC50 values for DCS as a function of arabinose induction level.

**GFP-ALR single cell induction.csv** \
Flow-cytometry data for the strain TESEC GFP-tagged Mtb ALR, collected at a range of arabinose concentrations. Used to show the induction profiles of the ALR target at single-cell resolution.

### Figure 2 A screen for targeted Mtb-ALR inhibitors identifies benazepril

**Prestwick Info 96 Format.csv** \
Names and properties of the compounds from the Prestwick library organized in 96-well format to match the screening data.

**TESEC ALR Screen DF 96.csv** \
Growth data for the TESEC ALR strain with the Prestwick library in 96-well format. Used to identify the initial set of hit compounds including benazepril. Intended to be human-readable and can also be imported as a Python DataFrame.

**TESEC ALR Screen Dict 96.pickle** \
The same growth data as the above file but packaged as a Python Dictionary to simplify computational analysis.

### Figure 3 Chemical-genetic drug sensitivity is significantly altered by target overexpression

**Prestwick Info 384 Format.csv** \
Names and properties of the compounds from the Prestwick library organized in 384-well format to match the screening data.

**TESEC ALR Screen DF 384.csv** \
Growth data for the TESEC ALR strain with the Prestwick library in 96-well format. Used to characterize the performance of the screen across a range of arabinose induction levels. Intended to be human-readable and can also be imported as a Python DataFrame.

**TESEC ALR Screen Dict 384.pickle** \
The same growth data as the above file but packaged as a Python Dictionary to simplify computational analysis.

### Figure 4 Chemical-genetic profiles indicate distinct mechanisms of action

**TESEC ChemGen Profiles.csv** \
Growth data for the TESEC ALR strain varying the concentrations of arabinose and 3 inhibitors: DCS, benazepril and amlexanox. Used to show that chemical-genetic profiles can suggest distinct mechanisms of action.

**TESEC Dala Rescue.csv** \
Growth data for the TESEC ALR strain treated with DCS or benazepril and rescued with D-alanine. This rescue activity is consistent with the idea that benazepril, like DCS, acts directly on ALR.

**TESEC Efflux Effect.csv** \
Growth data for the TESEC ALR strain and the TESEC Mtb ALR Efflux+ strain, in which the TolC efflux system is restored. Used to show the effect of efflux on TESEC sensitivity to DCS and benazepril.

### Figure 5 Antibiotic activity of benazepril against M. smegmatis and M. tuberculosis

**Msmeg Dala Rescue.csv** \
Growth data for M. smegmatis treated with DCS or benazepril and rescued with D-alanine.

**Msmeg IC50 Curves.csv** \
Growth Data for M. smegmatis treated with varying concentrations of DCS, benazepril, benazeprilat and DMSO. Used to confirm the activity of benazepril against myobacteria.

### Figure 6 Benazepril inhibits Mtb ALR in vitro.

**In Vitro Hanes Woolf.csv** \
Kinetic rates of purified ALR enzymes varying the concentration of both substrate (D-alaine) and inhibitor (benazepril). Used to differentiate competitive and noncompetitive mechanisms of inhibition.

**In Vitro Michaelis Menten.csv** \
Kinetic rates of purified ALR enzymes varying the concentration inhibitors DCS and benazepril. Used to find the best-fit inhibition constant under standard conditions.

### Style files
The .mplstyle files are used with matplotlib to set figure sizes and other properties.

**Kelly Colors.json** \
A set of hex color codes adapted from Kelly's 22 colors of maximum contrast (Kelly, K. Color Eng. 3(6), 1965)

### Python files

**requirements.txt** \
A list of python packages and versions used for these scripts. Created with pipreqs 0.4.10.

**pyenv.pex** \
A virtual environment supplying the required packages. Created with PEX 2.1.28.

## License

[Creative Commons Attribution-NonCommercial 4.0 International Public License](https://creativecommons.org/licenses/by-nc/4.0/)

## Acknowledgments

Thanks to the Bettencourt Schueller Foundation long term partnership, this work was partly supported by the CRI Research Fellowship to Edwin Wintermute. Additional support was provided by MSDAVENIR.
