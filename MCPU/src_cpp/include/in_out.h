#ifndef IN_OUT_H
#define IN_OUT_H

#include <stdio.h>

#include "define.h"
#include "globals.h"

void read_pdb_backbone(
    struct Simulation *sim,
    char pdb_name[], char res_name[5][4], double r_n[5][3], double r_a[5][3],
                       double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res,
                       int end_res);
void write_pdb_backbone(
    struct Simulation *sim,
    char pdb_name[], char res_name[5][4],
    Eigen::Matrix<double, 3, 5>& r_n,
    Eigen::Matrix<double, 3, 5>& r_a,
    Eigen::Matrix<double, 3, 5>& r_c,
    Eigen::Matrix<double, 3, 5>& r_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int stt_res, int end_res);

#endif