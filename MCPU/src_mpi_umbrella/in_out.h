#ifndef IN_OUT_H
#define IN_OUT_H

void read_pdb_backbone(char pdb_name[200], char res_name[5][4], double r_n[5][3], double r_a[5][3],
                       double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res,
                       int end_res);
void write_pdb_backbone(char pdb_name[200], char res_name[5][4], double r_n[5][3], double r_a[5][3],
                        double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res,
                        int end_res);

FILE *fpdb;

#endif