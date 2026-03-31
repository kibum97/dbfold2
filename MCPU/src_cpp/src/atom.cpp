#include "atom.h"

void CopyAtom(struct atom protein1, struct atom *protein2) {
    (*protein2).X            = protein1.X;
    (*protein2).Y            = protein1.Y;
    (*protein2).Z            = protein1.Z;
    (*protein2).xyz_int.x    = protein1.xyz_int.x;
    (*protein2).xyz_int.y    = protein1.xyz_int.y;
    (*protein2).xyz_int.z    = protein1.xyz_int.z;
    (*protein2).xyz.x        = protein1.xyz.x;
    (*protein2).xyz.y        = protein1.xyz.y;
    (*protein2).xyz.z        = protein1.xyz.z;
    (*protein2).res_num      = protein1.res_num;
    (*protein2).atomtype     = protein1.atomtype;
    (*protein2).smogtype     = protein1.smogtype;
    (*protein2).is_core      = protein1.is_core;
    (*protein2).is_designed  = protein1.is_designed;
    (*protein2).is_sidechain = protein1.is_sidechain;
    strcpy((*protein2).res, protein1.res);
    strcpy((*protein2).atomname, protein1.atomname);
    (*protein2).matrix        = protein1.matrix;
    (*protein2).sec_structure = protein1.sec_structure;
    return;
}