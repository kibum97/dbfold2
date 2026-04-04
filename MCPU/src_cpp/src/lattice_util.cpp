#include "lattice_util.h"

void CopyLatticeCoordinates(struct atom A, struct atom *B) {
    (*B).X         = A.X;
    (*B).Y         = A.Y;
    (*B).Z         = A.Z;
    (*B).xyz_int   = A.xyz_int;
    (*B).matrix    = A.matrix;
    return;
}

void UpdateLattice(
    struct Context *ctx, struct MCIntegrator *integrator,
    short rotate_natoms, short *rotate_atom) {
    int j;

    for (j = 0; j < rotate_natoms; j++) {
        temp_atom    = &ctx->native[rotate_atom[j]];
        temp_atom->X = PerBound(integrator, (int)(floor(temp_atom->xyz[0] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
        temp_atom->Y = PerBound(integrator, (int)(floor(temp_atom->xyz[1] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
        temp_atom->Z = PerBound(integrator, (int)(floor(temp_atom->xyz[2] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
        temp_atom->matrix    = &(ctx->the_matrix[temp_atom->X][temp_atom->Y][temp_atom->Z]);
        temp_atom->xyz_int = temp_atom->xyz.cast<int>();
    }

    return;
}

/* computes which lattice cell the given atom falls into */
/* also, finds the integer representation of the coordinates */
void FindLatticeCoordinates(
    struct Context *ctx, struct MCIntegrator *integrator,
    struct atom *ATOM) {
    (*ATOM).X         = PerBound(integrator, (int)(floor((*ATOM).xyz[0] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
    (*ATOM).Y         = PerBound(integrator, (int)(floor((*ATOM).xyz[1] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
    (*ATOM).Z         = PerBound(integrator, (int)(floor((*ATOM).xyz[2] * integrator->LATTICE_SIZE)) + integrator->HALF_MATRIX_SIZE);
    (*ATOM).xyz_int = (*ATOM).xyz.cast<int>();
    (*ATOM).matrix    = &(ctx->the_matrix[(*ATOM).X][(*ATOM).Y][(*ATOM).Z]);

    return;
}

void InitializeMatrix(
    struct Context *ctx, struct MCIntegrator *integrator
) {
    /*The Matrix is a  3D array, each dimension is given by MATRIX_SIZE = 20 (# of AAs)*/
    short i, j, k, a, b, c;

    ctx->the_matrix = (struct cell ***)calloc(integrator->MATRIX_SIZE, sizeof(struct cell **));
    for (i = 0; i < integrator->MATRIX_SIZE; i++) {
        ctx->the_matrix[i] = (struct cell **)calloc(integrator->MATRIX_SIZE, sizeof(struct cell *));
        for (j = 0; j < integrator->MATRIX_SIZE; j++)
            ctx->the_matrix[i][j] = (struct cell *)calloc(integrator->MATRIX_SIZE, sizeof(struct cell));
    }
    for (i = 0; i < integrator->MATRIX_SIZE; i++)
        for (j = 0; j < integrator->MATRIX_SIZE; j++)
            for (k = 0; k < integrator->MATRIX_SIZE; k++) {
                ctx->the_matrix[i][j][k].natoms = 0;
                ctx->the_matrix[i][j][k].X      = i;
                ctx->the_matrix[i][j][k].Y      = j;
                ctx->the_matrix[i][j][k].Z      = k;
                for (a = -1; a < 2; a++)
                    for (b = -1; b < 2; b++)
                        for (c = -1; c < 2; c++) {
                            ctx->the_matrix[i][j][k].neighbors[(a + 1) * 9 + (b + 1) * 3 + (c + 1)] =
                                &(ctx->the_matrix[PerBound(integrator, i + a)][PerBound(integrator, j + b)][PerBound(integrator, k + c)]);
                        }
            }

    return;
}

unsigned char PerBound(
    struct MCIntegrator *integrator,
    signed char X) {
    signed char value;

    if (X < 0) {
        value = X + integrator->MATRIX_SIZE;
        while (value < 0)
            value += integrator->MATRIX_SIZE;
        return value;
    } else if (X >= integrator->MATRIX_SIZE) {
        value = X - integrator->MATRIX_SIZE;
        while (value >= integrator->MATRIX_SIZE)
            value -= integrator->MATRIX_SIZE;
        return value;
    } else
        return X;
}
