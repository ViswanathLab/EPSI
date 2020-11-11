Instructions
Step 1: Import "high resolution anatomical reference image" using file 'read_one_image_ref.m'. Run the code after giving the folder number. Save the output (RefImage) as .mat file.
Step 2: Import "EPSI spectroscopic data" using file 'Open_BrukerEPSI_2dseq.m'. Run the code after inputing the folder number. The code has to be run in increments. Run until line 36. Then depending on the file (Magnitude or Phase data) run either line 38 (Magnitude data) or line 40 & 41 (Phase data). After you have run the code twice, once for magnitute and once for phase, save the final output (MetImage) as .mat file.
Step 3: Now that you have the anatomical reference image and the spectroscopic EPSI data as matlab matrixes, you are ready to further process your data (spectra zero fill & apodize, calculate metabolic maps and evaluate SNR) using 'CodeMatlab_Final.m'. Also in this script you can plot spectra and metabolic maps.

Example data attached
folder 2 - HR anatomical data (step 1)
folder 8 - EPSI metabolic data (step 2)

Files included (Matlab R2018a)
Main files:
read_one_image_ref.m
Open_BrukerEPSI_2dseq.m
CodeMatlab_Final.m
Function file:
readBrukerHeader.m








