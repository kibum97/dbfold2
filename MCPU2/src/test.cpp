#include <Eigen/Dense>
#include <iostream>
#include <vector>

class SubMatrixView {
public:
    SubMatrixView(Eigen::MatrixXd& matrix, const std::vector<int>& rowIndices)
        : matrix_(matrix), rowIndices_(rowIndices) {}

    Eigen::Ref<Eigen::VectorXd> row(int i) {
        return Eigen::Map<Eigen::VectorXd>(matrix_.row(rowIndices_[i]).data());
    }

    int rows() const {
        return rowIndices_.size();
    }

    int cols() const {
        return matrix_.cols();
    }

private:
    Eigen::MatrixXd& matrix_;
    std::vector<int> rowIndices_;
};

int main() {
    // Create a main matrix (N*3)
    Eigen::MatrixXd mainMatrix(8, 3);
    mainMatrix << 1, 2, 3,
                  4, 5, 6,
                  7, 8, 9,
                  10, 11, 12,
                  13, 14, 15,
                  16, 17, 18,
                  19, 20, 21,
                  22, 23, 24;

    // Indices of the rows to include in the submatrix
    std::vector<int> rowIndices = {1, 5, 3, 7};

    // Create a submatrix view
    SubMatrixView subMatrix(mainMatrix, rowIndices);

    // Print the submatrix
    std::cout << "Submatrix:\n";
    for (int i = 0; i < subMatrix.rows(); ++i) {
        std::cout << subMatrix.row(i).transpose() << std::endl;
    }

    // Modify the main matrix
    mainMatrix(5, 1) = 42;

    // Print the modified main matrix
    std::cout << "Modified main matrix:\n" << mainMatrix << std::endl;
    std::cout << mainMatrix.row(5) << std::endl;

    // Print the submatrix to see the change
    std::cout << "Submatrix after modification:\n";
    for (int i = 0; i < subMatrix.rows(); ++i) {
        std::cout << subMatrix.row(i).transpose() << std::endl;
    }

    return 0;
}