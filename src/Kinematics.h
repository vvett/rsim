#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "MassProperties.h"

namespace rsim {

    // defines motion of a body
    class Kinematics {
        private:
            std::vector<double> time;
            double mass;
            Vector3d cg_loc; // in body frame 

        public:

            void update();

    };

}
