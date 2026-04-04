#include "jac_local.h"
#include "tripep_closure.h"

void loop_Jacobian(Mat3 r_n, 
                   Mat3 r_ca, 
                   Mat3 r_c, 
                   double *Jac) {
    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        std::array<Vec3, 6> axis;
        std::array<Vec3, 6> pivot;
        int n = 0;
        // Compute Aces and Pivots
        for (int i = 0; i < 3; i++) {
            // phi
            axis[n] = (r_ca.col(i) - r_n.col(i)).normalized();
            pivot[n] = r_ca.col(i);
            n++;
            // psi
            axis[n] = (r_c.col(i) - r_ca.col(i)).normalized();
            pivot[n] = r_c.col(i);
            n++;
        }
        const Vec3& r_ca3 = r_ca.col(2);
        const Vec3& r_cac3 = axis[5];

        // Determine which components to use for Jacobian
        int c1, c2;
        if (std::abs(r_cac3.z()) < 1.0e-10) {
            c1 = 0; c2 = 2; // Use x and z
        } else {
            c1 = 0; c2 = 1; // Use x and y
        }

        // Compute determinant and Jacobian
        Eigen::Matrix<double, 5, 5> J5;
        J5.setZero();

        for (n = 0; n < 4; n++) {
            J5.block<1, 3>(n, 0) = axis[n].cross(r_ca3 - pivot[n]).cast<double>().transpose();
        }
        for (n = 0; n < 4; n++) {
            Vec3 ang = axis[n].cross(r_cac3);
            J5(n, 3) = ang[c1];
            J5(n, 4) = ang[c2];
        }
        double det = J5.determinant();
        
        // Successfully computed Jacobian, break out of loop
        if (std::abs(det) > 1.0e-10) {
            *Jac = 1.0e0 / fabs(det);
            break;
        }

        // Jacobian is too close to singular, perturb the input and try again
        Eigen::Quaternion<double> ran_q(
            threefryrand() - 0.5, // w
            threefryrand() - 0.5, // x
            threefryrand() - 0.5, // y
            threefryrand() - 0.5  // z
        );
        ran_q.normalize();
        Eigen::Matrix3d ran_U = ran_q.toRotationMatrix();
        for (int i = 0; i < 3; i++) {
            r_n.col(i)  = ran_U * r_n.col(i);
            r_ca.col(i) = ran_U * r_ca.col(i);
            r_c.col(i)  = ran_U * r_c.col(i);
        }
    }
    *Jac = -1.0e0; // Indicate total failure
}