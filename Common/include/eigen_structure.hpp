
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "datatype_structure.hpp"

using MatrixType = Eigen::Matrix<su2double, Eigen::Dynamic, Eigen::Dynamic>;
using QRdecomposition = Eigen::ColPivHouseholderQR<MatrixType>;

using VectorType = Eigen::Matrix<su2double, Eigen::Dynamic, 1>;
