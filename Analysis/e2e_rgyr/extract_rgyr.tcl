# Author: Nicholas Yerolemou - Apr 2024

# Usage
# 1. Load your data into VMD
# 2. Update MOL, output_Path, output_Name and sel
# 3. Run the script with 'source extract_e2e.tcl'




# -----EDIT HERE ------ #
set MOL "Pn23F_6RU_V2"
set output_Path "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/${MOL}/Analysis/rgyr/"
set output_Name "${MOL}_0_to_1000ns_rgyr.txt"
set sel [atomselect top "noh and not type SOD and not resname G2P ARHM and not resid 1 30 and not name C6 O2 O3 O6 and not index 520"]



set output_location "${output_Path}${output_Name}"
set outFileWriter [open $output_location "w"]


set ntot [molinfo top get numframes]
# measure the radius of gyration and save to output 
for { set frame 0 } { $frame < $ntot } { incr frame } {
	molinfo top set frame $frame

	set rgyr [measure rgyr $sel]

	#write output
	puts $outFileWriter $frame\t$rgyr
	}

# close output 
close $outFileWriter
