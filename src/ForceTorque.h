#pragma once
#include "defs.h"
#include "Frames.h"

namespace rsim {

    // provides a standard interface for forces and moments
    class ForceTorque {

        private:

            // Frame that the Force is applied in
            Frame* force_frame;
        
        public:
            // gets torque about certain point
            Vector3d getTorque(Vector3d& ref_point);
            
            // gets force
            Vector3d getForce();

            // gets 6x1 force torque about certain point: TODO - suggest inline this
            Vector3d getForceTorque(Vector3d& ref_point);

            // updates force
            void update();

    };
}