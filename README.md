Peptide Landscape Antigenic Epitope Alignment Utility (PLAtEAU)

A Python script for identifying, aligning, and quantifying antigenic epitopes displayed by MHCII molecules.

Authors: Miguel Alvaro-Benito, Eliot Morrison, Esam Abualrous, Benno Kuropka, and Christian Freund
(FU Berlin)

#####################

*** PLAtEAU is now available as a free web tool! It can be accessed at https://plateau.bcp.fu-berlin.de/

#####################

Instructions

This script requires the following libraries: numpy, scipy, itertools, pyteomics

1. Perform identification of immunopeptidome mass spectrometric data using the MaxQuant software with the parameters of your choice.

2. (optional) Filter the evidence.txt file to your specifications. Examples include removing potential contaminants, peptides identified with low confidence, or background peptides from IP controls. 

3. Copy the .fasta file used in (1) to your working folder with the filtered evidence.txt file. 

4. In this folder, run the script:
 python plateau-1.0.py -exp="experiment name" -evidence="evidence.txt" -fasta="fasta_file.fasta"
    
5. The final output file is called "experiment name"_core_epitopes_final_renorm.txt
The "Core Epitope" column represents the conserved antigenic epitopes. Core epitopes marked with * represent potential frame shifts and should be investigated closely.
The "Whole Epitope" column represents the entire primary sequence covered by the identified peptides. 
