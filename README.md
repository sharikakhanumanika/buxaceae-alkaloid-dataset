buxaceae-alkaloid-dataset

Computational dataset and analysis files for the investigation of natural alkaloids from the Buxaceae plant family as multi-target anti-Alzheimer's agents: A computational study

Research Overview

The study combines molecular docking, ADMET prediction, molecular dynamics (MD) simulations, and binding free-energy analysis to investigate the interaction and stability of selected Buxaceae alkaloids with AChE and BChE.

Contents
1. Molecular Docking
Docking analysis performed using PyRx (AutoDock Vina, blind docking)
Raw docking output files and binding affinity data are available from the corresponding author upon reasonable request — see Data Availability Note below
2. Molecular Dynamics Simulations

The repository contains processed MD simulation data and analysis results, including:

RMSD
RMSF
Radius of gyration (Rg)
Solvent accessible surface area (SASA)
Hydrogen-bond analysis
Principal component analysis (PCA)
Dynamic cross-correlation analysis (DCCM)
Free Energy Landscape (FEL)
Molecular Surface Area (MolSA)
Polar Surface Area (PolSA)
3. Binding Free-Energy Analysis
MM-GBSA results for selected protein–ligand complexes
4. Figures

Figures generated from the computational analyses are provided for visualization and interpretation of the results.

5. Scripts

Python and other analysis scripts used for processing and visualization of the computational data are included where applicable.

Software and Tools
Docking and Screening
AutoDock Vina (within PyRx) — molecular docking and virtual screening; blind docking with an exhaustiveness value of 8 and one docking run per protein–ligand pair
BIOVIA Discovery Studio and PyMOL — protein–ligand interaction visualization and interaction diagram generation
ADMET and Toxicity Prediction
SwissADME — prediction of blood–brain barrier (BBB) permeability using the BOILED-Egg model, gastrointestinal (GI) absorption, CYP inhibition, lipophilicity (LogP), water solubility (LogS), and Lipinski's Rule of Five
ProTox-3.0 — in silico acute toxicity prediction
Structure Preparation and Topology Generation
UCSF ChimeraX (v1.11.1) — protein structure preparation using the DockPrep utility
GROMACS pdb2gmx — protein topology generation
sort_mol2_bonds.pl utility — correction and preparation of ligand MOL2 bond information
SwissParam server — generation of CHARMM-compatible ligand topology and parameters
Bioactivity Prediction
PASS (Prediction of Activity Spectra for Substances; Way2Drug) — prediction of biological activity spectra for the selected compounds
Molecular Dynamics Simulation
GROMACS 2025.3 — molecular dynamics simulations using the CHARMM36 force field and TIP3P water model
Molecular Dynamics Trajectory Analysis
GROMACS trjconv, gmx hbond, gmx covar, and gmx anaeig modules — trajectory processing, hydrogen-bond analysis, covariance analysis, and principal component analysis (PCA)
MDAnalysis (Python library) — custom trajectory analysis and dynamic cross-correlation matrix (DCCM) calculation
Python (Matplotlib) — in-house script for three-dimensional free-energy landscape (FEL) visualization
Repository Contents
Figure/ — figures generated from the computational analyses
MD_analysis/ — RMSD, RMSF, Rg, SASA, hydrogen-bond, PCA, DCCM, MolSA, PolSA, and FEL analysis result files
MMGBSA/ — MM-GBSA binding free-energy results
scripts/ — Python scripts used for data analysis and visualization

Raw MD simulation trajectories (GROMACS input/output files) and raw molecular docking output are not hosted in this repository due to file size — see Data Availability Note below.

Molecular Dynamics Simulation Details
Simulation engine: GROMACS 2025.3
Force field: CHARMM36
Water model: TIP3P
Simulation length: 100 ns per trajectory
Number of independent trajectories: 42 (36 holo and 6 apo; triplicate simulations)
Time step: 2 fs
Temperature: 300 K
Pressure: 1 bar
Analysis
Root mean square deviation (RMSD)
Root mean square fluctuation (RMSF)
Radius of gyration (Rg)
Solvent accessible surface area (SASA)
Hydrogen-bond analysis
Principal component analysis (PCA)
Dynamic cross-correlation matrix (DCCM)
Free-energy landscape (FEL)
MM-GBSA binding free-energy analysis
Molecular surface area (MolSA)
Polar surface area (PolSA)
Data Availability Note

Raw molecular docking output files (PDBQT, logs) and raw MD simulation trajectories (GROMACS input/output files) are available from the corresponding author upon reasonable request, due to file size and format constraints. ADMET (SwissADME, ProTox-3.0) and PASS prediction outputs can be regenerated using the SMILES codes provided in the manuscript via the respective publicly available web servers.

Reproducibility

The repository provides processed computational data, analysis outputs, figures, and relevant scripts used in the study. Software versions and computational parameters are documented where applicable to facilitate reproducibility.

Citation

If you use the data or scripts from this repository, please cite: [Authors]. [Dataset title]. Zenodo. DOI: [DOI]
