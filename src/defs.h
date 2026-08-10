#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace rsim {

    using Matrix3d = Eigen::Matrix<double, 3, 3>;
    using Vector3d = Eigen::Matrix<double, 3, 1>;

    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

}