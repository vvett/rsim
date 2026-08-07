#include "MassProperties.h"
#include "ForceTorque.h"
#include "defs.h"
#include <Eigen/Geometry>

namespace rsim {

    class Kinematics {

        private:
            struct Frame {
                Vector3d position;
                Eigen::Quaterniond orientation;
            };


        public:
            double timestep;
            Vector3d acceleration;
            Vector3d velocity;
            Vector3d 


        
    };

    class {};

    class {};


}