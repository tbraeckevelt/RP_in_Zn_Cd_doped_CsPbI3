# RP_in_Zn_Cd_doped_CsPbI3
A repository for the paper about the creation of RP phases due to Zn and Cd doping in CsPbI3

This repository consists of three folders:

1. `struc_generation/`
    - Contains python scripts and initial xyz files for structure generation.
    - This folder contains three subfolders:
        - `perovskite_subst_scripts/`: Using a 4 by 4 by 4 supercell of the alpha phase, Zn and Cd doped perovskite are created by substiting Pb with Zn or Cd, or by substituting Cs with Zn+I or Cd+I.
        - `scripts_for_RP_subst/`: This folder contains 5 subfolders 
            ```
            - `RP_subst_1/`: Contains substitution scripts for RP phase 1.
            - `RP_subst_2/`: Contains substitution scripts for RP phase 2.
            - `RP_subst_3/`: Contains substitution scripts for RP phase 3.
            - `RP_subst_4/`: Contains substitution scripts for RP phase 4.
            - `RP_subst_5/`: Contains substitution scripts for RP phase 5.
            ```
        - `interface_creation_sciripts/`: Description of subfolder3.

2. `MD_simulations_and_analysis/`
    - Contains files and scripts for molecular dynamics simulations and analysis.
    - sfg

3. `Train_macemodel/`
    - Contains files and scripts related to training the MACE model.
    - sfg


    > **Note:** This repository contains additional structures and Python scripts for molecular dynamics simulations which are not reported in the paper. Additional structures include doped perovskite and an alternative one shift 110 RP phase (110_1shifta, 110_1shiftb is reported in the paper due to its better metastability).