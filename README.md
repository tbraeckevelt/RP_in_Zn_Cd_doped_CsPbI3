# RP_in_Zn_Cd_doped_CsPbI3
A repository for the paper about the creation of RP phases due to Zn and Cd doping in CsPbI3

This repository consists of three folders:

1. `struc_generation/`
    - Contains python scripts and initial xyz files for structure generation.
    - This folder contains three subfolders:
        - `perovskite_subst_scripts/`: Using a 4 by 4 by 4 supercell of the alpha phase, Zn and Cd doped perovskite are created by substiting Pb with Zn or Cd, or by substituting Cs with Zn+I or Cd+I.
        - `scripts_for_RP_subst/`: This folder contains 5 subfolders, each subfolder contains Zn and Cd doped RP structures which are created by substiting Pb with Zn or Cd, or by substituting 2Cs with Zn or Cd in a pure 336-atoms supercell of the corresponding RP structure.
            ```
            - `001_1shift/`: Scripts for RP phase with CsI layer along a 001 plane of the alpha phase. The perovskite layer are shifted along one direction of the plane.
            - `001_2shifts/`: Scripts for RP phase with CsI layer along a 001 plane of the alpha phase. The perovskite layer are shifted along the two directions of the plane.
            - `110_1shifta/`: Scripts for RP phase with CsI layer along a 110 plane of the alpha phase. The perovskite layer are shifted along direction a of the plane. This 'a' direction is a 001 vector of the alpha phase
            - `110_1shiftb/`: Scripts for RP phase with CsI layer along a 110 plane of the alpha phase. The perovskite layer are shifted along direction b of the plane. This 'b' direction is a 1-10 vector of the alpha phase. This is the one-shift RP phase structure for the Cd doped samples that are reported in the paper.
            - `110_2shifts/`: Scripts for RP phase with CsI layer along a 110 plane of the alpha phase. The perovskite layer are shifted along the two directions of the plane.
            ```
        - `interface_creation_sciripts/`: This folder contains the script Create_perovskite_RP_interface.py, which constructs interface structures between the perovskite and RP phase. You can define the interface plane, the CsI plane of the RP phase and the various thicknesses in each direction. The perovskite-RP interface can be interpreteded as two peroskite structures shifted with a Pb-I distance in all three directions, in between those two perovskite structures is an S-shaped CsI layer. To create a one-shift RP phase the shift between those two perovskite structures is only performed in two directions. However, in most cases this will lead to overlapping atoms of different elements.

2. `MD_simulations_and_analysis/`
    - Contains files and scripts for molecular dynamics simulations and analysis.
    - sfg

3. `Train_macemodel/`
    - Contains files and scripts related to training the MACE model.
    - sfg


    > **Note:** This repository contains additional structures and Python scripts for molecular dynamics simulations which are not reported in the paper. Additional structures include doped perovskite and an alternative one shift 110 RP phase (110_1shifta, 110_1shiftb is reported in the paper due to its better metastability).