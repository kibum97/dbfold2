# -*- tcl -*-
package require namdenergy
set seqsta STAXXX
set seqend ENDXXX
set tar_outfile "ENERGYOUTFILE"

mol new INPSFXXX
mol addfile INDCDXXX type dcd first 4 last 4 step 1

set outfile [open $tar_outfile w]
set sel0 [atomselect top "protein and noh"]

animate write pdb MIN_OUT_PDB_TEMP.pdb sel $sel0 mol 

mol delete all 

quit
