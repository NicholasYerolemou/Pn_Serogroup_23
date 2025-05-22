README
Author: Nicholas Yerolemou

This is an introduction to using the extract_Dihedrals_all.tcl script
The goal of this script is to automatically extact all the dihedrals of every glycosidic linkage in your molecule
based on details supplied by you about the first repeat unit in the molecule.

# Step 1 - Specifiy RESID's and Atom Names for each linkage type:

    **1.1** In the **edit here** section update/create selections 1->N.

    Each selection specifies a unique glycosidic linkage in your first RU.
    If you have 5 unique linkage types per RU then you should have 5 selections
    The resid's and atom names that make up the linkage are specified in the format:

    set sel1 {first residue RESID {"atom 1 name" "atom 2 name"} second residue RESID {"atom 3 name" "atom 4 name" "atom 5 name"} "linkage name for output file"}

    Here is an example of a bDGal(1->4)bLRha linkage

        set sel1 {2 {"H1" "C1"} 1 {"O4" "C4" "H4"} "bDGal_14_bLRha"}

        This is a linkage from bDGal (resid 2) To bLRha (resid 1) defined by atoms H1 C1 O4 C4 H4.
        Where H1 and C1 belong to the Gal residue and O4 C4 and H4 belong to the Rha residue.

        Note: The resid of bDGal (2) is higher than bLRha (1) even though the linkage starts at bDGal and ends at bLRha.
        This is due to the order residues are built when the pdb files are created.

        NB: You need to correctly specify which atoms belong to which residue.


    **1.2** Once you have created a selection for each unique linkage update the **selections** list.
    If you have 5 linkages the selections list should look as follows **set selections [list $sel1 $sel2 $sel3 $sel4 $sel5]**.

    NB: The sel that defines the glycosidic linkage that links repeating units together must be placed last in selections.
    There are only n-1 occurances of this linkage compared to n occurances of the others.

    In the case that you want to define "omega" and "epsilon" dihedrals.

    Simply define 6 or 7 atom names instead of 5

    Omega is defined as atoms (3 4 5 6) in the list
    Epsilon is defined as atoms (4 5 6 7) in the list

    For 6 atoms Phi, Psi and Omega files are produced.
    For 7 atoms Phi, Psi, Omega and Epsilon files are produced.

    Here is an example of Omega and Epsilon being defined

    set sel2 {3 {"H2" "C2" "O2" "P1"} 2 {"O3" "C3" "H3"} "Gro2P_(O->3)_bDGal"}

    This defines 7 atom names and so angles are defined as:
    Phi: "H2" "C2" "O2" "P1"
    Psi: "C2" "O2" "P1" "O3"
    Omega: "O2" "P1" "O3" "C3"
    Epsilon: "P1" "O3" "C3" "H3"




    **1.3** Update the "MOL", "PATH", "repeat_units", "residues_per_unit" variables at the beginning of the edit here section


    ---------------
    For Pn23B 6RU, it should look as follows:
    set repeat_units 6
    set residues_per_unit 4
    set sel1 {2 {"H1" "C1"} 1 {"O4" "C4" "H4"} "bDGal_14_bLRha"}
    set sel2 {3 {"H2" "C2" "O2" "P1"} 2 {"O3" "C3" "H3"} "Gro2P_(O->3)_bDGal"}
    set sel3 {4 {"H1" "C1"} 2 {"O4" "C4" "H4"} "bDGlc_14_bDGal"}
    set sel4 {5 {"H1" "C1"} 4 {"O4" "C4" "H4"} "bLRha_14_bDGlc"}


    set selections [list $sel1 $sel2 $sel3 $sel4]
    ---------------

# Step 2 - Validate your selections

    **2.1** Find the MAIN-START HERE section at the bottom of the script. uncomment step 1 and comment out step 2
    **2.2** Load your molecule with its dcd data into VMD and open the tcl console
    **2.3** Cd to the directory that holds this script and run it using **source extract_Dihedrals_All.tcl**


    Data will not yet be extracted but the createRepresentations function will be run.
    The createRepresentations function is there to help validate your selections from step 1 before extacting the dihedral data.
    It will create representations in red in vmd to show the linkages you have defined.
    **2.4** View your molecule in VMD and check that red representaions have been defined for every glycosidic linkage in your molecule. Check that the correct atoms for each linkage are defined.

    If this is correct contine, otherwise repeat step 1.
    Common mistakes: - Incorrect RESID order. The RESID's should be listed in the order they define the linkage, i.e. from RESID 1 to RESID 2 - Incorrect atom names for that residue. The names in the list after the RESID belong to that residue.
    Make sure you have not mixed up which atoms belong to which residue/RESID.

    If the linkages are all correct proceed to Step 3.

# Step 3 - Extract the data

    In the "MAIN-START HERE" section at the bottom of the script comment out Step 1 and uncomment Step 2.
    Run the script to extract the dihedral data.

    **Output**:
    The output files are written to the location defined by the "PATH" varibale.
    A file is created for each sel defined in step 1 and named acording to the filename given in that definition.
    The output is in the format **frame,Phi,Psi,Omega,Epsilon**. This is written as a header at the start of the file.
    If there are no values for Omega or epsilon the comma's will still be present just with no value between them e.g. "5,18.7,-12.9,,"

    If the molecule has multiple RU's there will be multuple occurances of the same linkage type, labelled A->Z.
    The data for occurance A is written first followed by occurance B.
    Before the data for each occurance a header is written.
        # Linkage Occurance C
        #PHI Atoms:177 176 171 169
        #PSI Atoms:176 171 169 170
    This tells you which occurance this is and what the index numbers are of the atoms that make up the Phi, Psi and potentially Omega, Epislon angles.

# Step 4 - Plot the data

    The plot_Dihedral_and_PMF.py script was made to parse the output of this script and plot it.

# Change Log

V2.4

- Fixed bug when creating representaions of dihedrals

V2.3

- All occurances of the same linkage type are written to the same file with headers in-between

V2.2 Change log

- All data Phi, Psi, Epsilon, Omega data for a linkage is written to a single file in format (frame,phi,psi,epsilon,omega)

V2.1 Change log

- prints atom index's for each dihedral at start of each file
- frameStart variable added to control start frame of loopy part of the name
- \_A is added to the first file if not alread

V2.0

- Input style changed to resid {atoms} resid {atoms}
- printOutResidueDetails function removed
- Linkage information printed at the start

# Notes:

#We cannot combine the selections in the calculateDihedralAngles function as getIndex/list does not preserve the order of the atom index's but sorts them in ascening order.
#Using seperate selections like we are forced to greatly increases run time
#set phi_selection [atomselect top "(resid $first_resid and name [lrange $atom_names 2 3]) || (resid $second_resid and name [lrange $atom_names 0 1])"]
#set psi_selection [atomselect top "(resid $first_resid and name [lrange $atom_names 2 end]) || (resid $second_resid and name [lrange $atom_names 1 1])"]
