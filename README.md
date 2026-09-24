# JLAB RG-C

This repository contains the code used for my study of Kaon SIDIS for single and double spin asymmetries at CLAS12 with the RGC experiment. 
Asymmetries observed: $A_{UL}^{\sin\phi}$, $A_{UL}^{\sin 2\phi}$, $A_{LL}$, $A_{LL}^{\cos\phi}$, $A_{LU}^{\sin\phi}$

## MAIN CODE

- **`rgc_analysis_1D.cpp`** 
Is the main code used to analyze the root file from the rgc campaign, which are converted and analyzed on the ifarm with my source code ('rgc.cpp'). This code is used to study single and double spin asymmetry in two dimension, on of $x_B-Q^2$ and the second of $z-P_{hT}$. The same analysis in 2+2D is inside 'rgc_analysis.cpp', but has to be modified since is not updated with all the correction of this code. Its produce a root file and a .csv table to use easily the extracted results.

- **`plot_asymmetries.py` / `plot_asymmetries.cpp`** 
Use the csv table from the previous code to produce a comparison among the different run period of the RGC experiment (summer22, fall22 and spring23). The py file produce a png output while the cpp a root file.

- **`rgc_multiplicity.cpp`**   
Used to calculate the multiplicity, and take as imput file a series of root dataset converted by myself on the ifarm. 

- **`rgc_new_mc_analysis.cpp`**  
Perform similar studies of 'rgc_analysis_1D.cpp' but with the MC production of RGC available from the collaboration on the farm. It's purpose is to observe the performance of an injected asymmetry in the MC tree, and observe our resolution in extracting it.

- **`rgc_toy_systematics.cpp`** 
Used to study the systematics by a toy model which repeat the asymmetry extractiuon, through a MLE, N times and observe the bias in the asymmetry produced by each step of reconstruction. At the moment is focused on four systematic: Tracking + Acceptance, PID, Purity & Contamination, and Bin Migration.

- **`analyze_toys.py`**  
Analyze the output from the toy code and produce different png plot and a .csv table with all the information.

- **`plot_asymmetries_wSys.py`**  
Produce a plot of the systematics overlapped with the actual asymmetries extracted by the real RGC run (from 'rgc_analysis_1D.cpp'), to highlights the impact of the systematics with the actual measurement.

- **`rich_rgc_analysis.cpp`**  
Was used to perform some studies of the RICH performance.