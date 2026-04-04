#include "getrms.h"
#include <Eigen/SVD>

float getrms(struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn) {
    int NFRAG = algn.NFRAG;
    
    // Count total aligned atoms to pre-allocate
    int n = 0;
    for (int i = 1; i <= NFRAG; i++) {
        n += (algn.seqptr[i].x2 - algn.seqptr[i].x1 + 1);
    }
    if (n == 0) return 0.0f;

    // Extract coordinates into 3xN Matrices (SVD is fastest on contiguous data)
    Eigen::Matrix3Xd P(3, n), Q(3, n);
    int col = 0;
    for (int i = 1; i <= NFRAG; i++) {
        for (int j = algn.seqptr[i].x1; j <= algn.seqptr[i].x2; j++) {
            int k = j - algn.seqptr[i].x1 + algn.structptr[i].x1;
            P.col(col) = struct1[j].CA;
            Q.col(col) = struct2[k].CA;
            col++;
        }
    }

    // Vectorized Centroid Subtraction (No loops!)
    P.colwise() -= P.rowwise().mean();
    Q.colwise() -= Q.rowwise().mean();

    // Calculate Covariance Matrix H
    Eigen::Matrix3d H = P * Q.transpose();

    // Run SVD
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    // Handle reflection (det(V*U^T) must be 1 for a rotation)
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    Eigen::Matrix3d eye = Eigen::Matrix3d::Identity();
    if (d < 0) eye(2, 2) = -1.0;

    Eigen::Matrix3d Rotation = svd.matrixV() * eye * svd.matrixU().transpose();

    // Calculate RMSD directly from SVD result
    // This is mathematically equivalent to the lowest eigenvalue method
    double ssd = P.squaredNorm() + Q.squaredNorm() - 2.0 * (svd.singularValues().array() * eye.diagonal().array()).sum();
    
    return std::sqrt(std::max(0.0, ssd / n));
}

float getrms_fast(struct backbone* struct1, struct backbone* struct2, 
                  struct alignment& algn, Eigen::Matrix3Xd& P, Eigen::Matrix3Xd& Q) {
    int n = 0;
    // 1. Efficiently fill only the columns needed for this specific alignment
    for (int i = 1; i <= algn.NFRAG; i++) {
        for (int j = algn.seqptr[i].x1; j <= algn.seqptr[i].x2; j++) {
            int k = j - algn.seqptr[i].x1 + algn.structptr[i].x1;
            // No allocation here, just copying 3 doubles into a pre-existing slot
            P.col(n) = struct1[j].CA;
            Q.col(n) = struct2[k].CA;
            n++;
        }
    }

    // 2. Use "TopLeftCorner" to only work on the active atoms
    // This allows you to use a large buffer for proteins of different sizes
    auto P_active = P.topLeftCorner(3, n);
    auto Q_active = Q.topLeftCorner(3, n);

    // 3. Vectorized Centroid Subtraction
    // rowwise().mean() is highly optimized in Eigen
    P_active.colwise() -= P_active.rowwise().mean();
    Q_active.colwise() -= Q_active.rowwise().mean();

    // 4. Compute Covariance H (Matrix multiplication is the fastest part of Eigen)
    Eigen::Matrix3d H = P_active * Q_active.transpose();

    // 5. Pre-allocate the SVD object to avoid internal allocations
    // Use JacobiSVD for 3x3 because it is extremely stable
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    // 6. Calculate Square Sum of Distances (SSD) using Singular Values
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    double s3 = (d < 0) ? -svd.singularValues()(2) : svd.singularValues()(2);
    
    double ssd = P_active.squaredNorm() + Q_active.squaredNorm() - 
                 2.0 * (svd.singularValues()(0) + svd.singularValues()(1) + s3);

    return std::sqrt(std::max(0.0, ssd / n));
}