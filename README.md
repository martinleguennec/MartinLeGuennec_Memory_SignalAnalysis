# MartinLeGuennec_Memory_SignalAnalysis
 
 All the explanations contained in this git as well as the code will be written in English because it is the convention for coding. The supplementary material, however, is written in French because it is part of the memory.
 This git contains 1 pdf and 2 repository, they will be detailed in this file.
 
 # Materiel_Supplementaire.pdf

This pdf document contains supplementary informations that couldn't fit in the memory.


# Matlab_EMG_signalTreatment

This repository contains the matlab code used to treat the EMG data. 

The EMG are not in this repository due to confidentiality issues. They are property of EuroMov DHM and can't be in public access.


The repository is structured as follows: 
- a PRG folder containing the matlab codes
- a RES folder containing the files produced by the matlab code

To run the code, simply run the main.m script.
All scripts are written according to the following guidelines: https://www.ee.columbia.edu/~marios/matlab/MatlabStyle1p5.pdf


# R_MVC prediction 

This folder contains the R script that builds the models used in the brief. It uses data from the matlab code presented above.

The folder has the following structure:
- A DAT folder containing the data file used by the code 
- A LE_GUENNEC.MARTIN.Rmd file containing the code executed to analyze the data with the models. The pdf and html files with the same name are created from the code, enabling you to obtain the results without having to run the code. 
- A R_MVCprediction.Rproj file, used to launch the Rmd script.
- A apa.csl and references.bib file, useful for creating the bibliography in the Rmd file