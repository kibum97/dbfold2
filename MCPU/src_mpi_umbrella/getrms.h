#ifndef GETRMS_H
#define GETRMS_H

#define ROTATE(a, i, j, k, l)        \
    g       = a[i][j];               \
    h       = a[k][l];               \
    a[i][j] = g - s * (h + g * tau); \
    a[k][l] = h + s * (g - h * tau);

struct fragments {
    int x1;
    int x2;
};

struct alignment {
    struct fragments seqptr[MAXFRAG];
    struct fragments structptr[MAXFRAG];
    int              NFRAG;
};

struct backbone {
    struct vector N;
    struct vector CA;
    struct vector C;
    struct vector CB;
    struct vector O;
};

float getrms(struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE],
             struct alignment algn);
int   kearsley(float a[5][5], struct backbone struct1[MAXSEQUENCE],
               struct backbone struct2[MAXSEQUENCE], struct alignment algn);
void  jacobi(float a[5][5], int n, float d[5], float v[5][5], int *nrot);

#endif