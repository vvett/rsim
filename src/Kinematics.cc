#include "Kinematics.h"

namespace rsim {

    // class that stores objects and methods to calculate rigid body kinematics

    /*
     * Calculate acceleration and velocity in body frame, then convert to an inertial frame
     * Propogate update to frame via RK4 integrator
     */


    Kinematics::Kinematics (Frame& frameBody, ForceTorque& forceTorque, RigidBody& rigidBody)
    : bodyFrame(frameBody), bodyForceTorque(forceTorque), bodyMassProperties(rigidBody) {};

    Vector6d Kinematics::solveAccels () {
        /*
         * Returns accels as column vector with body origin, body frame translational acceleration
         * first followed by body frame anguar acceleration. Solve via Newton Euler.
         */

        Vector6d accelVector;

        double mass = bodyMassProperties.mass();
        Vector3d cg_loc = bodyMassProperties.cg();
        Matrix3d moi_cg = bodyMassProperties.momentOfInertia();

        Vector6d ft_body = bodyForceTorque.getForceTorque(cg_loc);

        Matrix3d I3 = Matrix3d::Identity();

        Matrix3d cskew;
        cskew << 0, -cg_loc(2, 0), cg_loc(1, 0),
                 cg_loc(2, 0), 0, -cg_loc(0, 0),
                 cg_loc(1, 0), cg_loc(0, 0), 0;

        Matrix3d wskew;
        wskew << 0, -bodyAngularVelocity(2, 0), bodyAngularVelocity(1, 0),
                 bodyAngularVelocity(2, 0), 0, -bodyAngularVelocity(0, 0),
                 bodyAngularVelocity(1, 0), bodyAngularVelocity(0, 0), 0;
        
        Matrix6d inertia_tensor;
        Matrix3d IQ1 = mass * I3;
        Matrix3d IQ2 = -mass * cskew;
        Matrix3d IQ3 = mass * cskew; 
        Matrix3d IQ4 = moi_cg - mass * cskew * cskew;
        inertia_tensor << IQ1, IQ2,
                          IQ3, IQ4;

        Vector6d coriolis_vector;
        Vector3d CQ1 = mass * wskew * wskew * cg_loc;
        Vector3d CQ2 = wskew * (moi_cg - mass * cskew * cskew) * bodyAngularVelocity;

        // TODO Solve accels


        
        return accelVector;
    };
}