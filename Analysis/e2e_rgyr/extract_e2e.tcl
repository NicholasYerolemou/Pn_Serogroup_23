# Author: Nicholas Yerolemou - Apr 2024

# Usage
# 1. Load your data into VMD
# 2. Update MOL, output_Path, output_Name, atom_1_index, and atom_2_index
# 3. Run the script with 'source extract_e2e.tcl'

# -----EDIT HERE ------ #
set MOL "Pn23F_6RU_V2"
set output_Path "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/${MOL}/Analysis/e2e/"
set output_Name "${MOL}_0_to_1000ns_e2e.txt"
set atom_1_index 22
set atom_2_index 518



set output_location "${output_Path}${output_Name}"
set outFileWriter [open $output_location "w"]

set ntot [molinfo top get numframes]
# Measure the bond length and save to output 
for { set frame 0 } { $frame < $ntot } { incr frame } {
    molinfo top set frame $frame

    # ATOM INDEX NUMBERS FOR E2E & MOLID -----
	set bond [measure bond [list $atom_1_index $atom_2_index] molid 0]

    # Write output
    puts $outFileWriter $frame\t$bond
}

# Close output
close $outFileWriter
