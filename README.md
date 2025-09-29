# NIDER project
 Data and analyses for a collaborative meta-analysis on the Neurobiology of Individual Differences in Emotion Regulation.


# Reproducibility Notes
This project uses Matlab and R. All analysis scripts can be found in the "scripts" folder. Matlab analyses require installing the CanlabCore Toolbox, which is available on github. Its installation sectin details the required dependencies (e.g., SPM). Further, the "functions" folder of this repository has to be on the matlab path for custom functions. 
Packages for R analyses are tracked with the renv package (via the respective renv folder on the top level). Open the "NIDER_project.Rproj" file on the top level. Then, install the "renv" package and use "renv::restore()". Then, open and use the R scripts from within this Rstudio environment.
The folder "testAnalyses" contains tests on the meta-analysis functions.
The OS was Microsoft Windows 11 Home. Analyses were performed using Matlab v2023b, the CanlabCore Toolbox (last tested on commit 884ef4c), and R version 4.3.2. Installation time is about an hour, mostly depending on the matlab installation time.
