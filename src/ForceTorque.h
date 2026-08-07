#pragma once
#include <array>
#include <vector>

// provides a standard interface for forces and moments
class ForceTorque {
    
    public:
        // gets torque about certain point

        
        // gets force
        std::array<double, 3> getForce();

        // updates force
        void update();

};