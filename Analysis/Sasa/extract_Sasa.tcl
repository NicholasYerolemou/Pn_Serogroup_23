# Author: Nicholas Yerolemou

# Usage
# 1. Load your data into VMD
# 2. Update MOL, Path, probe, output_Name and hydrophobic_selection
# 3. Run the script with 'source extract_e2e.tcl'


#number of equilibration frames to skip
set skip 0
set molid 0
set probe 2.5

#Small probe = 1
#Medium probe = 1.4
#Large probe = 2.5


# -----EDIT HERE ------ #
set MOL "Pn23A_9RU"
set PATH "/home/nicholas-yerolemou/Documents/UCT/PhD/Simulation/Pn23/9RU/${MOL}/Analysis/Sasa/"
set output_Name "${MOL}_SASA_Large_Gro2P.txt" 

set hydrophobic_selection "resid 24"

#The name of every *hydrophobic* atom in the system
# set hydrophobic_selection "name H1 H2 H3 H4 H5 H32 H31 H51 H52 H61 H62 HT1 HT2 HT3 H63 HM HM3 HM2 HM1 H11 H12 C1 C2 C3 C4 C5 C6 CM CM2 CT C CTB CAB"

#The name of every *hydrophilic* atom in the system
# set hydrophilic_selection "name HO1 HO2 HO3 HO4 HO6 OA O O1 O2 O3 O4 O5 O6 O61 O62 OA N HN HP2 P1 OP3 OP4 OP2 HO5 O1B O2B NB"


# --------------------- #


#Setup output file
set output "${PATH}${output_Name}"
set outFileWriter [open $output "w"]


# Obtain total number of frames 
set ntot [molinfo $molid get numframes] 

#Loop through each frame
for { set frame $skip } { $frame < $ntot } { incr frame } { 
	molinfo $molid set frame $frame
	
	set sel [atomselect $molid  "segname CARB"]
	#Measure total sasa region on molecule
	set tot [measure sasa $probe $sel ]

	# -----Hydrophobic ------ #
	set hydrophobic [atomselect $molid $hydrophobic_selection]
	set sasa_phobic [measure sasa $probe $sel -restrict $hydrophobic ]

	#sasa expressed as a % (hydrophibic region)/totalsasa
	set sasa [expr $sasa_phobic/$tot*100]
	# ---------------------- #


	# # -----Hydrophilic ------ #
	# set hydrophilic [atomselect $molid $hydrophilic_selection]
	# set sasa_philic [measure sasa $probe $sel -restrict $hydrophilic ]

	# #sasa expressed as a % (hydrophibic region)/totalsasa
	# set sasa [expr $sasa_philic/$tot*100]
	# # ---------------------- #


	puts $outFileWriter "$frame\t$sasa"

}

close $outFileWriter