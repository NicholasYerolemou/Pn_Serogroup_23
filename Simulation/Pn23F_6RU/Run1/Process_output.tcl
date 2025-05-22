#load the PDB and DCD files
animate read pdb ../Pn23F_6RU_V2_Min_H2O_Na.pdb beg 0 end -1 skip 1 waitfor all
animate read dcd Pn23F_6RU_V2_Min_H2O_Na_run_1.dcd beg 0 end -1 skip 100 waitfor all

#Select everything except water
set x [atomselect top "not water"] ; 

#Output new PDB and DCD files
animate write dcd Pn23F_6RU_0_to_200ns.dcd beg 0 end -1 sel $x
