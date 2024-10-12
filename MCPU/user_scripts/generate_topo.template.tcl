# -*- tcl -*-
package require psfgen
topology top_all27_prot_lipid.inp
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

mol new INPDBXXX

set sel0 [atomselect top protein]
$sel0 writepdb temp.pdb
mol delete all

segment PROT {pdb temp.pdb}
coordpdb temp.pdb PROT
guesscoord
writepdb OUTPDBXXX
writepsf OUTPSFXXX

exit
