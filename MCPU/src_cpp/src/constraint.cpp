#include "constraint.h"

#include <fstream>

void Init_frag_rmsd_constraints(
    struct Simulation *sim
) {
    char             str[1000];
    int              frag_idx;
    int              frag_x1;
    int              frag_x2;
    struct fragments constraint_seqptr[MAXFRAG];
    struct fragments constraint_structptr[MAXFRAG];

    std::ifstream file(sim->rmsd_constraint_file);
    int   nfrag = 0;
    if (!file.is_open()) {
        fprintf(sim->STATUS, "Could not open file %s \n", sim->rmsd_constraint_file);
        exit(1);
    } else {
        fprintf(sim->STATUS, "Successfully opened constraint file %s \n", sim->rmsd_constraint_file);
        fflush(sim->STATUS);
    }
    while (file.getline(str, 1000)) {
        if (str[0] != '#') {
            if (sscanf(str, "%d %d %d", &frag_idx, &frag_x1, &frag_x2) != 3)
                continue;
            constraint_seqptr[frag_idx + 1].x1    = frag_x1;
            constraint_structptr[frag_idx + 1].x1 = frag_x1;
            constraint_seqptr[frag_idx + 1].x2    = frag_x2;
            constraint_structptr[frag_idx + 1].x2 = frag_x2;
            nfrag++;
        }
    }
    file.close();

    constraint_align.NFRAG = nfrag;
    for (int i = 1; i <= nfrag; i++) {
        constraint_align.seqptr[i]    = constraint_seqptr[i];
        constraint_align.structptr[i] = constraint_structptr[i];
    }

    return;
}

float Compute_frag_rmsd_constraint_energy(
    struct System *sys,
    struct backbone  struct1[MAXSEQUENCE],
    struct backbone  struct2[MAXSEQUENCE],
    struct alignment algn) {
    // This function will compute the energy of the fragment rmsd constraint
    float frag_rmsd = getrms(struct1, struct2, algn);
    return sys->k_constraint * frag_rmsd * frag_rmsd;
}

int Argmin(
    struct System *sys,
    float distances[]) {
    int argmin = 0;
    int i;
    for (i = 0; i < sys->N_constraints; i++) {
        if (distances[i] < distances[argmin]) {
            argmin = i;
        }
    }
    return argmin;
}

void Read_constraints(
    struct System *sys,
    struct Simulation *sim,
    struct Context *ctx
) {
    sys->constraint_array = (int **)calloc(
        sys->Max_N_constraints,
        sizeof(
            int *));  // Soo...we are allocating an array of pointers ( Max_N_constraints of them)
    // The calloc function is saying that we want to allocate memory for an array of
    // Max_N_constraints integer pointers The (int **) at the beginning is saying that malloc should
    // return a pointer to an array of pointers to ints (the memory for this array has been
    // allocated)

    sys->constraint_weights = (float *)calloc(sys->Max_N_constraints,
                                         sizeof(float));  // Return a pointer to an array of floats

    // char* filename = "Sample_constraint.txt";

    char  str[1000];
    int   aa;
    int   bb;
    float cc;
    int   iii = 0;
    std::ifstream file(sim->constraint_file);
    if (!file.is_open()) {
        fprintf(sim->STATUS, "Could not open file %s \n", sim->constraint_file);
        exit(1);
    } else {
        fprintf(sim->STATUS, "Successfully opened constraint file %s \n", sim->constraint_file);
        fflush(sim->STATUS);
    }
    while (file.getline(str, 1000)) {
        if (str[0] != '#') {
            // fprintf(STATUS, "%s", str);
            // fflush(STATUS);
            if (sscanf(str, "%d %d %f", &aa, &bb, &cc) != 3)
                continue;
            // fprintf(STATUS, "Did sscanf \n");
            // fflush(STATUS);
            // fprintf(STATUS, "The a value is %d while the iii value is %i \n", aa, iii);
            // fflush(STATUS);

            sys->constraint_array[iii] = (int *)calloc(
                2, sizeof(int));  // This is saying that the ith element of constrain_array is a
                                  // pointer to memory for 2 integers

            sys->constraint_array[iii][0] = aa;
            sys->constraint_array[iii][1] = bb;
            sys->constraint_weights[iii]  = cc;
            iii++;
            // fprintf(STATUS, "Hi this is line %i \n", iii);
            // fflush(STATUS);
        }
    }

    file.close();
    sys->N_constraints = iii;
    if (sys->N_constraints > 0) {
        fprintf(sim->STATUS, "We have %i constraints. The 1st constraint is %i, %i with weight %f \n",
                sys->N_constraints, sys->constraint_array[0][0], sys->constraint_array[0][1],
                sys->constraint_weights[0]);
        fflush(sim->STATUS);
    }

    sys->distances =
        (float *)calloc(sys->N_constraints, sizeof(float));  // Will contain distance between i and j
    ctx->disulfide_pairs_attempt = (int *)calloc(
        sys->N_constraints,
        sizeof(int));  // Note that when you allocate an array with calloc, the default value for
                       // the elements is zero. This is in contrast to malloc

    return;
}

float Compute_constraint_energy(
    struct Context *ctx, struct System *sys,
    const struct residue *residues, const struct atom *atoms) {
    float constraint_energy = 0;
    int   i, res1, res2, atomi, atomj;
    float weight, xi, yi, zi, xj, yj, zj, dist;
    float sum_of_weights = 0;

    int   min_ind;
    float min_dis;

    // disulfides = (int*)calloc(N_cysteines,sizeof(int)); //Is residue i participating in
    // disulfide?

    /* int disulfides[1000] = {0}; //Tels you whether each residue is participating in a disulfide
     */

    ctx->mean_constraint_distance = 0;

    for (i = 0; i < sys->N_constraints; i++) {
        res1   = sys->constraint_array[i][0];
        res2   = sys->constraint_array[i][1];
        weight = sys->constraint_weights[i];

        atomi = residues[res1].CB;
        if (atomi == -999) {
            /* GLY CB gets assigned index -999.
             *  In this case we will use the CA.
             */
            atomi = residues[res1].CA;
        }
        
        atomj = residues[res2].CB;
        if (atomj == -999) {
            /* GLY CB gets assigned index -999.
             *  In this case we will use the CA.
             */
            atomj = residues[res2].CA;
        }
        atomi = residues[res1].CA;
        atomj = residues[res2].CA;
        dist = (atoms[atomi].xyz - atoms[atomj].xyz).norm();


        /* Note we don't use this var */
        sys->distances[i] = dist;
        /* Note we don't use this var */
        ctx->disulfide_pairs_attempt[i] = 0;  // for now, reset the value of disulfide..we will figure
                                         // out later if a disulfide is present

        weight = sys->constraint_weights[i];
        /* flat-bottom harmonic which returns -k_constraint (* weight)
         * for distances within the region. */
        if (dist <= 2.9) {
            constraint_energy += weight * sys->k_constraint * ((dist - 2.9) * (dist - 2.9) - 1);
        } else if (dist <= 8.) {
            constraint_energy += weight * sys->k_constraint * -1;
        } else { /* (dist > 8.) */
            constraint_energy += weight * sys->k_constraint * ((dist - 8.) * (dist - 8.) - 1);
        }
    }

    /* if (mcstep % MC_PRINT_STEPS == 0) */
    /* { */
    /* 	fprintf(STATUS, "actual disulfide %d %d %d %f\n", */
    /* 		min_ind, I, J, min_dis); */
    /* } */

    /* distances[min_ind] = INFINITY; */
    /* min_ind = Argmin(distances); */
    /* min_dis = distances[min_ind]; */

    // if (mcstep==0){
    // fprintf(STATUS, "Let's check on our variable disulfide_pairs_attempt");
    // for (zzz=0; zzz<N_constraints; zzz++){
    //   if (disulfide_pairs_attempt[zzz]==1){
    //   	fprintf(STATUS, "%i and %i \n", constraint_array[zzz][0],constraint_array[zzz][1]);
    //   }
    //   }
    // }
    return constraint_energy;
}

///////////////////////////////////////////////

float Compute_constraint_energy_old(
    struct Context *ctx, struct System *sys,
    const struct residue *residues, const struct atom *atoms) {
    float constraint_energy = 0;
    int   i, res1, res2, atomi, atomj;
    float weight, xi, yi, zi, xj, yj, zj, dist;
    float sum_of_weights = 0;

    ctx->mean_constraint_distance = 0;

    for (i = 0; i < sys->N_constraints; i++) {
        res1   = sys->constraint_array[i][0];
        res2   = sys->constraint_array[i][1];
        weight = sys->constraint_weights[i];

        // To do, check that neither i nor j are already forming a disulfide bond before proceeding

        atomi = residues[res1].CA;
        atomj = residues[res2].CA;
        dist = (atoms[atomi].xyz - atoms[atomj].xyz).norm();

        constraint_energy = constraint_energy + weight * sys->k_constraint * (dist - 5) * (dist - 5);

        // Add here: if this dist is between 4 and 6, add residues i and j to the list of residues
        // that are already part of a disulfide bond
        sum_of_weights           = sum_of_weights + weight;
        ctx->mean_constraint_distance = ctx->mean_constraint_distance + weight * dist;

        // if (mcstep % MC_PRINT_STEPS == 0) {
        // fprintf(STATUS, "%Distance  is %5.3f \n", dist );
        // fflush(STATUS);
        // }
    }
    ctx->mean_constraint_distance = ctx->mean_constraint_distance / sum_of_weights;

    //	if (mcstep % MC_PRINT_STEPS == 0) {
    // fprintf(STATUS, "Constraint energy is %5.3f \n", constraint_energy);
    // fflush(STATUS);
    //}

    return constraint_energy;
}
