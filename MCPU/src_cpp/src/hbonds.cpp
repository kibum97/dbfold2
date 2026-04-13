#include "hbonds.h"
#include <string.h>
#include "init.h"
#include "pdb_util.h"
#include "tripep_closure.h"
#include "vector.h"

#include <vector>
#include <unordered_map>
#include <string>

float HydrogenBonds(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
) {
    int   i, j, res_n, res_o;
    float e;

    e = 0;
    for (res_n = 1; res_n < top->nresidues - 1; res_n++)
        for (res_o = 1; res_o < top->nresidues - 1; res_o++) {
            i = ctx->native_residue[res_n].N;
            j = ctx->native_residue[res_o].O;
            if (ctx->data[i][j].check_hbond) {
                ctx->data[i][j].hbond = ctx->data[j][i].hbond = CheckHBond(ctx, sys, top, i, j);
                if (d_memory < HB_INNER * HB_INNER)
                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 0;
                else if ((d_memory > HB_INNER * HB_INNER) && (d_memory < HB_CUTOFF * HB_CUTOFF))
                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 1;
                else if (d_memory > HB_CUTOFF * HB_CUTOFF)
                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 3;
            }
            if (ctx->data[i][j].check_hbond)
                if (ctx->data[i][j].hbond != sys->NO_HBOND) {
                    if (abs(ctx->native[i].res_num - ctx->native[j].res_num) > 4) {
                        if (d_memory < HB_INNER * HB_INNER) {
                            e += beta_favor *
                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                       [GetAminoNumber(ctx->native[j].res)] *
                                 ctx->hbond_E[ctx->data[i][j].hbond];
                        } else
                            e += HB_PENALTY * beta_favor *
                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                       [GetAminoNumber(ctx->native[j].res)] *
                                 ctx->hbond_E[ctx->data[i][j].hbond];
                    } else {
                        if (d_memory < HB_INNER * HB_INNER)
                            e += sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                       [GetAminoNumber(ctx->native[j].res)] *
                                 ctx->hbond_E[ctx->data[i][j].hbond];
                        else
                            e += HB_PENALTY *
                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                       [GetAminoNumber(ctx->native[j].res)] *
                                 ctx->hbond_E[ctx->data[i][j].hbond];
                    }
                }
        }
    e = e / 1000.0 * RDTHREE_CON;
    return e;
}

float FoldHydrogenBonds(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
) {
    int   i, j, k, res_n;
    int   d;
    float e;

    e = 0;

    for (res_n = 1; res_n < top->nresidues - 1; res_n++) {
        if (ctx->native_residue[res_n].res == "LNK") {
            continue;
        }
        i               = ctx->native_residue[res_n].N;
        temp_atom       = &ctx->native[i];
        temp_xyz_int    = temp_atom->xyz_int;
        temp_cell3      = ctx->prev_native[i].matrix;
        temp_cell_array = temp_atom->matrix->neighbors;
        for (d = 0; d < 27; d++) {
            temp_cell = temp_cell_array[d];
            if (temp_cell->natoms) {
                k = 0;
                O = temp_cell->natoms;
                A = temp_cell->atom_list;
                while (k < O) { /* this is O not zero */
                    j = A[k];

                    if (ctx->data[i][j].check_hbond) {
                        ctx->data[i][j].hbond = ctx->data[j][i].hbond = CheckHBond(ctx, sys, top, i, j);
                        if (d_memory < HB_INNER * HB_INNER)
                            ctx->data[i][j].closehb = ctx->data[j][i].closehb = 0;
                        else if ((d_memory > HB_INNER * HB_INNER) &&
                                 (d_memory < HB_CUTOFF * HB_CUTOFF))
                            ctx->data[i][j].closehb = ctx->data[j][i].closehb = 1;
                        else if (d_memory > HB_CUTOFF * HB_CUTOFF)
                            ctx->data[i][j].closehb = ctx->data[j][i].closehb = 3;

                        if (ctx->data[i][j].hbond != sys->NO_HBOND) {
                            if (abs(ctx->native[i].res_num - ctx->native[j].res_num) > 4) {
                                if (d_memory < HB_INNER * HB_INNER) {
                                    e += beta_favor *
                                         sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                               [GetAminoNumber(ctx->native[j].res)] *
                                         ctx->hbond_E[ctx->data[i][j].hbond];
                                } else
                                    e += HB_PENALTY * beta_favor *
                                         sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                               [GetAminoNumber(ctx->native[j].res)] *
                                         ctx->hbond_E[ctx->data[i][j].hbond];
                            } else {
                                if (d_memory < HB_INNER * HB_INNER) {
                                    e += sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                               [GetAminoNumber(ctx->native[j].res)] *
                                         ctx->hbond_E[ctx->data[i][j].hbond];
                                } else
                                    e += HB_PENALTY *
                                         sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                               [GetAminoNumber(ctx->native[j].res)] *
                                         ctx->hbond_E[ctx->data[i][j].hbond];
                            }
                        }
                    }
                    k++;
                }
            }
        }

        if (temp_cell_array[13] != temp_cell3->neighbors[13]) {
            for (d = 0; d < 27; d++) {
                k = 0;
                while (k < 27 && (temp_cell_array[k] != temp_cell3->neighbors[d]))
                    k++;
                if (k == 27) {  // if an old neighbor cell is no longer a neighbor cell...
                    temp_cell = temp_cell3->neighbors[d];
                    if (temp_cell->natoms) {
                        k = 0;
                        O = temp_cell->natoms;
                        A = temp_cell->atom_list;
                        while (k < O) {
                            j = A[k];

                            if (ctx->data[i][j].check_hbond) {
                                ctx->data[i][j].hbond = ctx->data[j][i].hbond = CheckHBond(ctx, sys, top, i, j);
                                if (d_memory < HB_INNER * HB_INNER)
                                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 0;
                                else if ((d_memory > HB_INNER * HB_INNER) &&
                                         (d_memory < HB_CUTOFF * HB_CUTOFF))
                                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 1;
                                else if (d_memory > HB_CUTOFF * HB_CUTOFF)
                                    ctx->data[i][j].closehb = ctx->data[j][i].closehb = 3;

                                if (ctx->data[i][j].hbond != sys->NO_HBOND) {
                                    if (abs(ctx->native[i].res_num - ctx->native[j].res_num) > 4) {
                                        if (d_memory < HB_INNER * HB_INNER) {
                                            e += beta_favor *
                                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                                       [GetAminoNumber(ctx->native[j].res)] *
                                                 ctx->hbond_E[ctx->data[i][j].hbond];
                                        } else
                                            e += HB_PENALTY * beta_favor *
                                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                                       [GetAminoNumber(ctx->native[j].res)] *
                                                 ctx->hbond_E[ctx->data[i][j].hbond];
                                    } else {
                                        if (d_memory < HB_INNER * HB_INNER) {
                                            e += sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                                       [GetAminoNumber(ctx->native[j].res)] *
                                                 ctx->hbond_E[ctx->data[i][j].hbond];
                                        } else
                                            e += HB_PENALTY *
                                                 sys->seq_hb[ctx->helix_sheet][GetAminoNumber(ctx->native[i].res)]
                                                       [GetAminoNumber(ctx->native[j].res)] *
                                                 ctx->hbond_E[ctx->data[i][j].hbond];
                                    }
                                }
                            }

                            k++;
                        }
                    }
                }
            }
        }
    }
    e = e / 1000.0 * RDTHREE_CON;
    return e;
}

void InitializeHydrogenBonding(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    struct Simulation *sim
) {
    int   i, j, k, l, match, skip;
    int   i7, i6, i5, i4, i3, i2, i1;
    int   jj1, jj2, jj3, jj4, jj5, jj6, jj7;
    int   D, D2, D3, D4, D5, D6, D7, A, A2, A3, A4, A5, A6, A7, A8;
    int   E;
    FILE *DATA;
    FILE *fseq_hb;
    char  line[250];
    char  line_seq_hb[250], file_seq_hb[250];
    char  res[4], res2[4], donor[4], donor2[4], donor3[4], donor4[4], donor5[4], donor6[4],
        donor7[4];
    char acceptor[4], acceptor2[4], acceptor3[4], acceptor4[4], acceptor5[4], acceptor6[4],
        acceptor7[4], acceptor8[4];
    int   offset[15];
    int   tmp_type, tmp_a, tmp_b;
    float tmp_val;


    float min_D, max_D, D_int;

    // if ((ctx->hbonds = (struct hydrogen_bond **)calloc(top->natoms, sizeof(struct hydrogen_bond *))) ==
    //     NULL) {
    //     fprintf(sim->STATUS, "ERROR: Failed to generate hbonds\n");
    //     exit(1);
    // }
    // for (i = 0; i < top->natoms; i++)
    //     ctx->hbonds[i] = (struct hydrogen_bond *)calloc(top->natoms, sizeof(struct hydrogen_bond));
    ctx->hbonds.clear();
    ctx->hbonds.resize(top->natoms);

    strcpy(file_seq_hb, "/n/home01/kibumpark/pkg/dbfold2/MCPU/src_cpp/config_files/seq_dep_hb_mu_low.energy");
    std::cout << "Opening the file: " << file_seq_hb << std::endl;
    fprintf(sim->STATUS, "Opening the file: %s\n", file_seq_hb);
    if ((fseq_hb = fopen(file_seq_hb, "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", file_seq_hb);
        exit(1);
        fclose(fseq_hb);
    }
    while (fgets(line_seq_hb, 250, fseq_hb) != NULL) {
        sscanf(line_seq_hb, "%d %d %d %f %*s%*s", &tmp_type, &tmp_a, &tmp_b, &tmp_val);
        if (sys->SEQ_DEP_HB == 0) {
            sys->seq_hb[tmp_type][tmp_a][tmp_b] = 1.;
        } else {
            sys->seq_hb[tmp_type][tmp_a][tmp_b] = (-1.) * tmp_val;
        }
    }
    std::cout << "DEBUG: Finished reading seq_dep_hb_mu_low.energy, first entry: " << sys->seq_hb[0][0][0] << std::endl;
    if ((DATA = fopen(sim->hydrogen_bonding_data.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->hydrogen_bonding_data);
        exit(1);
    }
    std::cout << "Opening the file: " << sim->hydrogen_bonding_data << std::endl;
    const auto& config = DEFAULT_HBOND_TOPOLOGY;

    for (int i = 0; i < top->nresidues; i++) {
        for (int j = 0; j < top->nresidues; j++) {
            // Check the skip parameter
            if (std::abs(i - j) < config.skip) continue;

            // Check boundary conditions for donors
            bool valid_bounds = true;
            for (int off : config.donor_offsets) {
                if (i + off < 0 || i + off >= top->nresidues) { valid_bounds = false; break; }
            }
            if (!valid_bounds) continue;

            // Check boundary conditions for acceptors
            for (int off : config.acceptor_offsets) {
                if (j + off < 0 || j + off >= top->nresidues) { valid_bounds = false; break; }
            }
            if (!valid_bounds) continue;

            // Atom number lookups for donors and acceptors
            try {
                std::array<int, 7> D_atoms;
                std::array<int, 8> A_atoms;

                for (size_t k = 0; k < 7; ++k) {
                    int idx = MatchAtomname(config.donors[k]); 
                    if (idx < 0) throw std::runtime_error("Atom not found");
                    D_atoms[k] = ctx->native_residue[i + config.donor_offsets[k]].atomnumber[idx];
                }
                for (size_t k = 0; k < 8; ++k) {
                    int idx = MatchAtomname(config.acceptors[k]); 
                    if (idx < 0) throw std::runtime_error("Atom not found");
                    A_atoms[k] = ctx->native_residue[j + config.acceptor_offsets[k]].atomnumber[idx];
                }

                int D = D_atoms[0];
                int A = A_atoms[0];

                // Populate the sparse maps symmetrically
                for (auto [src, dst] : {std::pair{A, D}, std::pair{D, A}}) {
                    auto& bond = ctx->hbonds[src][dst];
                    
                    bond.max_D2 = config.max_D * config.max_D;
                    bond.min_D2 = config.min_D * config.min_D;
                    bond.D_int  = config.D_int;

                    // Directly copy the resolved pointers into the array
                    for (size_t k = 0; k < 7; ++k) bond.donors[k]   = &(ctx->native[D_atoms[k]]);
                    for (size_t k = 0; k < 8; ++k) bond.acceptors[k] = &(ctx->native[A_atoms[k]]);
                }

                ctx->data[A][D].check_hbond = 1;
                ctx->data[D][A].check_hbond = 1;

            } catch (const std::exception&) {
                // If an atom isn't found in the topology (MatchAtomname fails),
                // we silently skip this specific bond, matching the original C behavior.
            }
        }
    }
    std::cout << "DEBUG: Finished reading hydrogen bonding data file." << std::endl;

    //  n_hbond_int = (int)((max_D-min_D)/D_int)+1;
    i7 = 3 * LTOR * LTOR * LTOR * LTOR * LANG * LANG;
    i6 = LTOR * LTOR * LTOR * LTOR * LANG * LANG;
    i5 = LTOR * LTOR * LTOR * LANG * LANG;
    i4 = LTOR * LTOR * LANG * LANG;
    i3 = LTOR * LANG * LANG;
    i2 = LANG * LANG;
    i1 = LANG;
    if ((ctx->hbond_E = (short *)calloc(i7, sizeof(short))) == NULL) {
        fprintf(sim->STATUS, "ERROR: Failed to generate hbond_E\n");
        exit(1);
    }
    sys->NO_HBOND = -1;
    for (jj1 = 0; jj1 < 3; jj1++)
        for (jj2 = 0; jj2 < LTOR; jj2++)
            for (jj3 = 0; jj3 < LTOR; jj3++)
                for (jj4 = 0; jj4 < LTOR; jj4++)
                    for (jj5 = 0; jj5 < LTOR; jj5++)
                        for (jj6 = 0; jj6 < LANG; jj6++)
                            for (jj7 = 0; jj7 < LANG; jj7++)
                                ctx->hbond_E[jj1 * i6 + jj2 * i5 + jj3 * i4 + jj4 * i3 + jj5 * i2 +
                                        jj6 * i1 + jj7] = 0;

    std::cout << "DEBUG: Finished initializing hbond_E array, starting to read energy values from file." << std::endl;

    std::string file = sim->hbond_energy_file;
    if ((sim->DATA = fopen(file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", file.c_str());
        exit(1);
    }

    fprintf(sim->STATUS, "Opening the file: %s\n", file);
    while (fgets(line, 250, sim->DATA) != NULL) {
        sscanf(line, "%d %d %d %d %d %d %d %d %*s%*s", &jj1, &jj2, &jj3, &jj4, &jj5, &jj6, &jj7,
               &E);
        ctx->hbond_E[jj1 * i6 + jj2 * i5 + jj3 * i4 + jj4 * i3 + jj5 * i2 + jj6 * i1 + jj7] = E;
        //    fprintf(STATUS, "%2d %2d %2d %2d %7d %5d\n", jj1, jj2, jj3, jj4,
        //    jj1*i3+jj2*i2+jj3*i1+jj4, E);
        //  fprintf(STATUS, "%d %d %d %d %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n",
        //  r1*20*i4+r2*i4+jj1*i3+jj2*i2+jj3*i1+jj4,
        // E, r1, r2, d1, d2, d3, d4, as, i1, jj1, jj2, jj3, jj4, r1*20*i4, r2*i4, jj1*i3, jj2*i2,
        // jj3*i1, jj4);
    }
    fclose(sim->DATA);

    //  min_E=0;
    //  for (i=0; i<20*20*i4; i++)
    //    if (hbond_E[i]<min_E)
    //      min_E = hbond_E[i];
    //  for (i=0; i<20*20*i4; i++)
    //    hbond_E[i]*=(hydrogen_bond/min_E);
    //  fprintf(STATUS, "hey %d %d %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n", r1, r2, d1, d2,
    //  d3, d4, i1, jj1, jj2, jj3, jj4, r1*20*i4, r2*i4, jj1*i3, jj2*i2, jj3*i1, jj4);

    std::cout << "DEBUG: Finished reading hydrogen bonding energy values from file." << std::endl;

    return;
}

int MatchAtomname(std::string_view name) {
    // .find() is the safest and fastest way to look up in a map
    if (auto it = ATOM_NAME_MAP.find(name); it != ATOM_NAME_MAP.end()) {
        return it->second; // Return the explicitly linked integer
    }    
    return -999; // Not found
}

long int CheckHBond(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top,
    int A, int B) {
    struct HydrogenBond *cur_hbond;
    Vec3         V3, H;
    Vec3         tmp1, tmp2, tmp3, tmp4;
    Vec3         plane1, plane2;
    Vec3         bisect1, bisect2;
    float                 d_HO;
    float                 d_CA_n0, d_CA_0p, d_CA_np, d_CA_00;
    float                 ang_Dphi, ang_Dpsi, ang_Aphi, ang_Apsi, ang_PH, ang_bH;
    float                 ang_PCA, ang_bCA, ang_PCAant, ang_bCAant;
    float                 ang_CACA;
    float                 min1, min2, min3;

    int res_dif;
    //  int helix_sheet = 0;
    int      jj1, jj2, jj3, jj4, jj5, jj6, jj7;
    long int xxx;

    cur_hbond = &(ctx->hbonds[A][B]);

    res_dif = abs((cur_hbond->donors[0])->res_num - (cur_hbond->acceptors[0])->res_num);

    if (((cur_hbond->donors[0])->res_num == 0) || ((cur_hbond->acceptors[0])->res_num == 0))
        return sys->NO_HBOND;
    if (((cur_hbond->donors[0])->res_num == top->nresidues - 1) ||
        ((cur_hbond->acceptors[0])->res_num == top->nresidues - 1))
        return sys->NO_HBOND;

    // estimate position of hydrogen atom, H
    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[1])->xyz, &H);
    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[2])->xyz, &V3);
    Add(V3, &H);
    Normalize(&H);
    Inverse(&H);
    Add((cur_hbond->donors[0])->xyz, &H);

    // determine if h...O distance is less than cutoff for hbonds
    d_HO     = D2(H, (cur_hbond->acceptors[0])->xyz);
    d_memory = d_HO;
    if (d_HO > HB_CUTOFF * HB_CUTOFF)
        return sys->NO_HBOND;

    // distances between CA atoms, p = previous, 0 = same, n = next residue
    d_CA_n0 = D2((cur_hbond->donors[4])->xyz, (cur_hbond->acceptors[2])->xyz);
    d_CA_0p = D2((cur_hbond->donors[1])->xyz, (cur_hbond->acceptors[5])->xyz);
    d_CA_np = D2((cur_hbond->donors[4])->xyz, (cur_hbond->acceptors[5])->xyz);
    d_CA_00 = D2((cur_hbond->donors[1])->xyz, (cur_hbond->acceptors[2])->xyz);
    min1    = min(d_CA_n0, d_CA_np);
    min2    = min(d_CA_0p, d_CA_00);
    min3    = min(min1, min2);
    //  fprintf(STATUS, "%d %f %f %f %f %f\n", res_dif, sqrt(d_CA_n0), sqrt(d_CA_0p), sqrt(d_CA_np),
    //  sqrt(d_CA_00), sqrt(min(min1, min2)));
    if (res_dif == 4) {
        ctx->helix_sheet = 0;
        if ((min1 > 5.8 * 5.8) || (min2 > 5.8 * 5.8))
            return sys->NO_HBOND;
        if (min3 > 5.5 * 5.5)
            return sys->NO_HBOND;
    }
    if (res_dif > 4) {
        if ((min1 > 6.0 * 6.0) || (min2 > 6.0 * 6.0))
            return sys->NO_HBOND;
        if (min3 > 5.4 * 5.4)
            return sys->NO_HBOND;
        // if secondary structure is predicted as helix, don't form h-bond when res_dif>4
        if ((top->secstr[(cur_hbond->donors[0])->res_num] == 'H') ||
            (top->secstr[(cur_hbond->acceptors[0])->res_num] == 'H'))
            return sys->NO_HBOND;
    }

    MakeVector((cur_hbond->donors[2])->xyz, (cur_hbond->donors[0])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[1])->xyz, &tmp2);
    MakeVector((cur_hbond->donors[1])->xyz, (cur_hbond->donors[3])->xyz, &tmp3);
    ang_Dphi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    ang_Dphi *= rad2deg;
    MakeVector((cur_hbond->donors[5])->xyz, (cur_hbond->donors[4])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[2])->xyz, &tmp2);
    MakeVector((cur_hbond->donors[2])->xyz, (cur_hbond->donors[0])->xyz, &tmp3);
    ang_Dpsi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    ang_Dpsi *= rad2deg;

    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp1);
    MakeVector((cur_hbond->acceptors[3])->xyz, (cur_hbond->acceptors[5])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[5])->xyz, (cur_hbond->acceptors[6])->xyz, &tmp3);
    ang_Aphi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    ang_Aphi *= rad2deg;
    MakeVector((cur_hbond->acceptors[4])->xyz, (cur_hbond->acceptors[2])->xyz, &tmp1);
    MakeVector((cur_hbond->acceptors[2])->xyz, (cur_hbond->acceptors[1])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp3);
    ang_Apsi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    ang_Apsi *= rad2deg;
    ang_Dphi += 180.;
    ang_Dpsi += 180.;
    ang_Aphi += 180.;
    ang_Apsi += 180.;
    if (res_dif == 4)
    // if we are in an alpha-helical conformation, but beta-sheet is predicted in sec_str file, no
    // h-bond
    {
        if ((top->secstr[(cur_hbond->donors[0])->res_num] == 'E') ||
            (top->secstr[(cur_hbond->acceptors[0])->res_num] == 'E')) {
            if ((ang_Dphi < 180.) && (ang_Dpsi < 180.))
                return sys->NO_HBOND;
            if ((ang_Aphi < 180.) && (ang_Apsi < 180.))
                return sys->NO_HBOND;
        }
        if ((top->secstr[(cur_hbond->donors[0])->res_num] == 'L') ||
            (top->secstr[(cur_hbond->acceptors[0])->res_num] == 'L')) {
            if ((ang_Dphi < 180.) && (ang_Dpsi < 180.))
                return sys->NO_HBOND;
            if ((ang_Aphi < 180.) && (ang_Apsi < 180.))
                return sys->NO_HBOND;
        }
    }
    if (res_dif > 4) {
        // don't form long-range h-bonds if we occupy quadrants on the right of ramachandran plot
        if (ang_Dphi > 150.)
            return sys->NO_HBOND;
        if ((ang_Dpsi > 30.) && (ang_Dpsi < 210.))
            return sys->NO_HBOND;
        if (ang_Aphi > 150.)
            return sys->NO_HBOND;
        if ((ang_Apsi > 30.) && (ang_Apsi < 210.))
            return sys->NO_HBOND;
        // don't form long-range h-bonds if L conformation is predicted
        if (top->secstr[(cur_hbond->donors[0])->res_num] == 'L')
            return sys->NO_HBOND;
        if (top->secstr[(cur_hbond->acceptors[0])->res_num] == 'L')
            return sys->NO_HBOND;
    }

    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[6])->xyz, &tmp1);
    MakeVector((cur_hbond->acceptors[7])->xyz, (cur_hbond->acceptors[5])->xyz, &tmp3);
    Normalize(&tmp1);
    Normalize(&tmp3);
    ang_CACA = Angle(tmp1, tmp3);
    ang_CACA *= rad2deg;
    if (res_dif > 4) {
        if (ang_CACA < 90)
            ctx->helix_sheet = 1;
        else
            ctx->helix_sheet = 2;
    }

    MakeVector((cur_hbond->donors[1])->xyz, (cur_hbond->donors[0])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[1])->xyz, (cur_hbond->donors[3])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[2])->xyz, (cur_hbond->acceptors[4])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[2])->xyz, (cur_hbond->acceptors[1])->xyz, &tmp4);
    CrossProduct(tmp1, tmp2, &plane1);
    CrossProduct(tmp3, tmp4, &plane2);
    ang_PCA = Angle(plane1, plane2);
    ang_PCA *= rad2deg;

    MakeVector((cur_hbond->donors[1])->xyz, (cur_hbond->donors[0])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[1])->xyz, (cur_hbond->donors[3])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[2])->xyz, (cur_hbond->acceptors[4])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[2])->xyz, (cur_hbond->acceptors[1])->xyz, &tmp4);
    Normalize(&tmp1);
    Normalize(&tmp2);
    Normalize(&tmp3);
    Normalize(&tmp4);
    bisect(tmp1, tmp2, &bisect1);
    bisect(tmp3, tmp4, &bisect2);
    ang_bCA = Angle(bisect1, bisect2);
    ang_bCA *= rad2deg;

    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[5])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[2])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[5])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[5])->xyz, (cur_hbond->acceptors[6])->xyz, &tmp4);
    CrossProduct(tmp1, tmp2, &plane1);
    CrossProduct(tmp3, tmp4, &plane2);
    ang_PCAant = Angle(plane1, plane2);
    ang_PCAant *= rad2deg;

    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[5])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[4])->xyz, (cur_hbond->donors[2])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[5])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[5])->xyz, (cur_hbond->acceptors[6])->xyz, &tmp4);
    Normalize(&tmp1);
    Normalize(&tmp2);
    Normalize(&tmp3);
    Normalize(&tmp4);
    bisect(tmp1, tmp2, &bisect1);
    bisect(tmp3, tmp4, &bisect2);
    ang_bCAant = Angle(bisect1, bisect2);
    ang_bCAant *= rad2deg;
    ang_Dphi = ang_PCA;
    ang_Dpsi = ang_bCA;
    ang_Aphi = ang_PCAant;
    ang_Apsi = ang_bCAant;

    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[2])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[1])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[2])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp4);
    CrossProduct(tmp1, tmp2, &plane1);
    CrossProduct(tmp3, tmp4, &plane2);
    ang_PH = Angle(plane1, plane2);
    ang_PH *= rad2deg;

    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[2])->xyz, &tmp1);
    MakeVector((cur_hbond->donors[0])->xyz, (cur_hbond->donors[1])->xyz, &tmp2);
    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[2])->xyz, &tmp3);
    MakeVector((cur_hbond->acceptors[1])->xyz, (cur_hbond->acceptors[3])->xyz, &tmp4);
    Normalize(&tmp1);
    Normalize(&tmp2);
    Normalize(&tmp3);
    Normalize(&tmp4);
    bisect(tmp1, tmp2, &bisect1);
    bisect(tmp3, tmp4, &bisect2);
    ang_bH = Angle(bisect1, bisect2);
    ang_bH *= rad2deg;
    //   fprintf(STATUS, "%d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", helix_sheet, ang_Dphi,
    //   ang_Dpsi, ang_Aphi, ang_Apsi, ang_PH, ang_bH); fprintf(STATUS, "%d %6.1f %6.1f %6.1f %6.1f
    //   %6.1f %6.1f\n", helix_sheet, ang_Dphi-180, ang_Dpsi-180, ang_Aphi-180, ang_Apsi-180,
    //   ang_PH, ang_bH);

    jj1 = ctx->helix_sheet;
    jj2 = (int)(ang_Dphi / ang3_int);
    jj3 = (int)(ang_Dpsi / ang3_int);
    jj4 = (int)(ang_Aphi / ang3_int);
    jj5 = (int)(ang_Apsi / ang3_int);
    jj6 = (int)(ang_PH / ang_int);
    jj7 = (int)(ang_bH / ang_int);

    xxx = jj7 + jj6 * LANG + jj5 * LANG * LANG + jj4 * LTOR * LANG * LANG +
          jj3 * LTOR * LTOR * LANG * LANG + jj2 * LTOR * LTOR * LTOR * LANG * LANG +
          jj1 * LTOR * LTOR * LTOR * LTOR * LANG * LANG;
    //  fprintf(STATUS, "%ld\n", xxx);
    return jj7 + jj6 * LANG + jj5 * LANG * LANG + jj4 * LTOR * LANG * LANG +
           jj3 * LTOR * LTOR * LANG * LANG + jj2 * LTOR * LTOR * LTOR * LANG * LANG +
           jj1 * LTOR * LTOR * LTOR * LTOR * LANG * LANG;
}
