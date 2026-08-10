#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "MassProperties.h"
#include "Frames.h"
#include "ForceTorque.h"

namespace rsim {

    // class that stores objects and methods to calculate rigid body kinematics

    /*
     * Calculate acceleration and velocity in body frame, then convert to an inertial frame
     * Propogate update to frame via RK4 integrator
     */
    class Kinematics {

        private:

            // time history vector
            std::vector<double> time;

            // body frame
            Frame& bodyFrame;

            // associated ForceTorque object
            ForceTorque& bodyForceTorque;

            // associated MassProperties object
            RigidBody& bodyMassProperties;


            // body frame kinematic properties
            Vector3d bodyAngularVelocity = {0, 0, 0};
            Vector3d bodyAngularAcceleration = {0, 0, 0};
            Vector3d bodyVelocity;
            Vector3d bodyAcceleration;

            // solves for angular and translational accel in body frame
            Vector6d solveAccels();

        public:
            
            Kinematics (Frame& frameBody, ForceTorque& forceTorque, RigidBody& rigidBody);
            
            void update();
    };
}