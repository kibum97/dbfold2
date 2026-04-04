#include "rotate.h"

#include <math.h>



void DoRotation(
    struct Context *ctx,
    int a, int b, int c, int d, Float delta_angle, short rotate_natoms,
    short *rotate_atom) {
    Vec3 pivot, axis;
    pivot = ctx->native[b].xyz;
    axis = (ctx->native[c].xyz - pivot).normalized();
    Eigen::Matrix3d rot_mat = Eigen::AngleAxisd(delta_angle, axis).toRotationMatrix();
    for (int j = 0; j < rotate_natoms; j++) {
        auto& current_xyz = ctx->native[rotate_atom[j]].xyz;
        current_xyz = rot_mat * (current_xyz - pivot) + pivot;
    }
    return;
}