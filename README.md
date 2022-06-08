![project schematic](https://github.com/FALL-ML/materials-discovery/blob/main/project_schematic.png)

# materials-discovery
Semi-supervised learning code + CAVD + BVSE + QE-NEB for identifying prospective ion conductors.

<p><b>Overview: Semi-supervised Learning.</b> The semi-supervised learning pipeline is split into six sub steps which are documented sequentially in the following notebooks:</p>

1. (Step1_Data_Pull_and_Clean.inyb): Pulls data from Materials Project and the ICSD folder. Merges structures that are contained in both databases. Charge decoration is applied to structures for enhanced compbatility with the descriptor transformations. Structure simplifications are applied. 
2. (Step2_Labeling_Data.inyb): Pulls in the labels csv and merges them with the structures from step 1.
3. (Step3_Generating_Descriptors.inyb): Applies the descriptor transformations to the structures (or simplified structures). 
4. (Step4_Screening_Descriptors.inyb): Not all structures are compatible with all descriptor transformations. The descriptor space is screened to identify a subset of structures that is compatible with all descriptors. 
5. (Step5_Agglomerative_Clustering.inyb): Agglomerative clustering is applied to the descriptor space. The resulting sets of 2-300 clusters are labeled with the room temperature ionic conductivity data. The Ward variance is calculated for the conductivity labels. 
6. (Step6_Comparing_Models.inyb): Code to visualize and compare each agglomerative clustering model. 

<p><b>Overview: CAVD+BVSE+NEB.</b> To expedite Nudged Elastic Band (NEB) calculations in Quantum Espresso, initials pathways were guessed by using Crystal Analysis by Voronoi Decomposition (CAVD) and Bond Valence Site Energy (BVSE) calculations. The notebook working_cavd_bvse_neb.inyb documents the entire process, all the way through performing an NEB calculation in Quantum Espresso. Qualitatively, the following steps are performed:</p>

1. Check that an input structure is ordered. Otherwise order it. 
2. Determine suitable hyperparamters by running self-consistent field (SCF) calculations on the input structure.
    1. kmesh - the k-space mesh.
    1. Ecutwfc - the wavefunction kinetic energy cutoff.
    1. Ecutrho - the charge density and potential kinetic energy cutoff. 
3. With suitable hyperparameters determined, relax the input structure to ~zero pressure.
4. Run the voronoi decomposition algorithm (CAVD) to find void pathways through the structure. 
5. Run the bond valence site energy (BVSE) calculation to estimate the electrostatic force field. 
6. Merge the CAVD and BVSE results to identify low-activation energy pathways. 
7. Select a pathway to examine. 
8. Generate an initial and final image from the pathway, by removing the relevant mobile ion. 
9. Relax the initial and final image to ~zero pressure. 
10. Create intermediate images by combining the CAVD+BVSE results with interpolation:
    1. Framework ions are interpolated between the initial and final images.
    1. Mobile ions are placed along the CAVD+BVSE pathway. 
11. Run & interperet Nudged Elastic Band on the images. 

***
The CAVD and BVSE libraries were originally generated by the Siqi Shi group. The code herein is altered for compability with Quantum Espresso. For VASP compatibility, their repo can be found here: https://gitee.com/shuhebing/cavd/tree/release


***

<p><b>Jupyter-lab Notes:</b> This jupyter notebook makes use of multiple jupyterlab extensions:</p>

* @jupyterlab/toc - for table of contents management
* @aquirdturtle/collapsible_headings - for collapsing sections

<p>Two additional extensions are used for git management:</p>

* @jupyterlab/git - for github integration
* @jupyterlab/nbdime-jupyterlab - for viewing revision history
