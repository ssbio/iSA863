# iSA863
Genome-scale metabolic model of Staphylococcus aureus USA300_FPR3757 (iSA863)

This is an updated and comprehensive genome-scale metabolic model of the methicillin-resistant human pathogen S. aureus  USA300_FPR3757. This repository contains necessary input files and scripts to simulate Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), Essentiality Analysis, Growmatch and other relevant optimization tool on the model in General Algebraic Modeling System (GAMS). 

Steps to follow when running a GAMS code:

Step 1: Create a GAMS scripts for the desired algorithm. Information on the syntax, options and guidelines for GAMS can be found in the GAMS user guide at https://www.gams.com/latest/docs/UG_Tutorial.html 
Step 2: Make sure that the current working directory is the directory that contains the file with the GAMS code. 
Step 3: When it finishes running, GAMS will output the results in files (if asked for) and state the location of these files.
Step 4: Open the file ending in .lst to find any errors in the GAMS code. Errors are marked with a dollar sign ($) in the .lst file. Descriptions of GAMS error messages and how to resolve them can be found at https://www.gams.com/latest/docs/UG_FixingErrors.html. 
Step 5: Open the file ending in .txt with a text editor to view the results of the GAMS code.
