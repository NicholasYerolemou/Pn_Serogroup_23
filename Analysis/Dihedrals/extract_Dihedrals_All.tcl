# V2.4
#See README_Extract_Dihedrals_All for usage and implementation details


#------------EDIT HERE------------#
set MOL "Pn23bb"
set molid 0

#Path to where you want the files saved
set PATH "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/Dihedrals/200_to_1000ns/"


set repeat_units 6
set residues_per_unit 3
puts "Your settings are:"
puts "$repeat_units repeat units, $residues_per_unit residues per unit"

set sel1 {2 {"H1" "C1"} 1 {"O4" "C4" "H4"} "bDGal_14_bLRha"}
set sel2 {3 {"H1" "C1"} 2 {"O4" "C4" "H4"} "bDGlc_14_bDGal"}
set sel3 {4 {"H1" "C1"} 3 {"O4" "C4" "H4"} "bLRha_14_bDGlc"}



set selections [list $sel1 $sel2 $sel3]

#---------------------------------#


set num_linkages_per_unit [llength $selections]
set alphabet "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#if bDGal_14_bLRha is the first linkage, All_selections is a 2d array holding all the occurances of bDGal_14_bLRha at pos 0. 
set All_Selections {} 

#loop through each linkage type and create selections for every time that linkage occurs
for {set n 0} {$n < $num_linkages_per_unit} {incr n} {

    set sel [lindex $selections $n]

    #Holds all the occurances of this linkage type
    set temp_selections_list [list $sel]
    set first_resid [lindex $sel 0]
    set second_resid [lindex $sel 2]
    set first_residue_atoms [lindex $sel 1]
    set second_residue_atoms [lindex $sel 3]


    # Get the names of the two residues forming the linkage
    set res1_atoms [atomselect top "resid ${first_resid}"] 
    set res2_atoms [atomselect top "resid ${second_resid}"]

    set resName1 [lindex [$res1_atoms get resname] 0]
    set resName2 [lindex [$res2_atoms get resname] 0]

    puts "Linkage from $resName1 to $resName2 with atoms $first_residue_atoms $second_residue_atoms"


    #loop through each repeat unit
    for {set i 1} {$i < $repeat_units} {incr i} {
        
        #Skip linkage because we have reached the terminal residue
        if {$n == [expr {$num_linkages_per_unit - 1}] && $i == [expr {$repeat_units - 1}]} {
            continue
            }


        set first_resid [expr {$first_resid + $residues_per_unit}]
        set second_resid [expr {$second_resid + $residues_per_unit}]

        set suffix [string index $alphabet [expr {$i}]]
        set selTemp1 [list $first_resid $first_residue_atoms $second_resid $second_residue_atoms $suffix]

        lappend temp_selections_list $selTemp1
    }
    lappend All_Selections $temp_selections_list
}


#------------CREATE REPS TO CHECK SELECTION------------#
proc createRepresentations {All_Selections} {
    # Get the number of representations
    set numReps [molinfo top get numreps]
    # Delete each old representation
    for {set i 0} {$i < $numReps} {incr i} {
        mol delrep top top
    }

    # Loop through selections
    foreach selection $All_Selections {
        foreach sel $selection {
            # Extract residue numbers and linkage atoms
            set first_resid [lindex $sel 0]
            set second_resid [lindex $sel 2]
            
            set first_residue_atoms [lindex $sel 1]
            set second_residue_atoms [lindex $sel 3]

            #create selection for this linkage
            set selection_string "(resid $first_resid and name $first_residue_atoms) || (resid $second_resid and name $second_residue_atoms)"


            # Apply new representation based on selection
            mol selection $selection_string
            mol representation Licorice 0.3 12
            mol color colorID 1
            mol addrep top 
        }
    }

    
    #add rep of whole molecule
    mol selection "not water and not segname ION"
    mol representation Licorice 0.2 12
    mol color ResName
    mol addrep top

}


#------------CALCULATE AND WRITE DIHEDRALS------------#
proc calculateDihedralAngles {selections PATH MOL molid} {

    #Open file for writing
    set file_name [lindex [lindex $selections 0] 4]
    set output_file [open "${PATH}${MOL}_${file_name}_Dihedrals.txt" a]

    #Write header
    puts $output_file "#Frame,Phi,Psi,Omega,Epsilon"

    foreach sel $selections {

        # Get details for this linkage
        set first_resid [lindex $sel 0]
        set second_resid [lindex $sel 2]
        set first_residue_atoms [lindex $sel 1]
        set second_residue_atoms [lindex $sel 3]
        set subtype [lindex $sel 4]

        # Write a header for the current subtype
        puts $output_file "#Linkage Occurrence $subtype"

        # Total number of atoms forming this linkage
        set num_atoms [expr {[llength $first_residue_atoms] + [llength $second_residue_atoms]}]

        # Loop through each frame and calculate + write dihedrals
        set numFrames [molinfo top get numframes]
        set startFrame 0
        for {set frame $startFrame} {$frame < $numFrames} {incr frame} {
            # Set the frame
            molinfo top set frame $frame 

            # Selection for this linkage's dihedrals
            set atoms {}
            foreach name $first_residue_atoms {
                set tempAtomIndex [[atomselect top "resid $first_resid and name $name"] get index]
                lappend atoms $tempAtomIndex
            }
            foreach name $second_residue_atoms {
                set tempAtomIndex [[atomselect top "resid $second_resid and name $name"] get index]
                lappend atoms $tempAtomIndex
            }

            # Initialize variables for dihedrals
            set epsilon ""
            set omega ""

            # Calculate additional dihedrals if applicable
            if {$num_atoms == 7} {
                set omega [measure dihed "[lindex $atoms 2] [lindex $atoms 3] [lindex $atoms 4] [lindex $atoms 5]" molid $molid]
                set epsilon [measure dihed "[lindex $atoms 3] [lindex $atoms 4] [lindex $atoms 5] [lindex $atoms 6]" molid $molid]
            } elseif {$num_atoms == 6} {
                set omega [measure dihed "[lindex $atoms 2] [lindex $atoms 3] [lindex $atoms 4] [lindex $atoms 5]" molid $molid]
            }

            # Calculate Phi and Psi dihedrals
            set phi [measure dihed "[lindex $atoms 0] [lindex $atoms 1] [lindex $atoms 2] [lindex $atoms 3]" molid $molid]
            set psi [measure dihed "[lindex $atoms 1] [lindex $atoms 2] [lindex $atoms 3] [lindex $atoms 4]" molid $molid]

            #Write which atoms make up this linkage
            if {$frame==$startFrame} {

                puts $output_file "#PHI Atoms:[lindex $atoms 0] [lindex $atoms 1] [lindex $atoms 2] [lindex $atoms 3]"
                puts $output_file "#PSI Atoms:[lindex $atoms 1] [lindex $atoms 2] [lindex $atoms 3] [lindex $atoms 4]"
                
                if {$num_atoms>5} {
                    puts $output_file "#OMEGA Atoms:[lindex $atoms 2] [lindex $atoms 3] [lindex $atoms 4] [lindex $atoms 5]"
                }
                if {$num_atoms>6} {
                    puts $output_file "#EPSILON Atoms:[lindex $atoms 3] [lindex $atoms 4] [lindex $atoms 5] [lindex $atoms 6]"
                }
            }

            # Write the frame number and dihedrals to the file (use empty string if Epsilon/Omega are not calculated)
            puts $output_file "$frame,$phi,$psi,$omega,$epsilon"
        }
    }
    # Close the file
    close $output_file
}




#------------MAIN-START HERE------------#

#Step 1.
# createRepresentations "$All_Selections"

#Step 2.
# Loops through each linkage type and writes out the data for each occurance of that type
foreach selections $All_Selections {
    calculateDihedralAngles "$selections" "$PATH" "$MOL" "$molid"
}


