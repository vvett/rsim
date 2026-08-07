#pragma once
#include "ForceTorque.h"
#include <vector>
#include <array>
#include <string>
#include "defs.h"

namespace rsim {
    
    class Vehicle {
        private:

            // array listing possible vehicle states
            std::vector<std::string> _vehicle_states = {"rail", "thrusting", "coasting", "descending", "drogue_deploy", "main_deploy", "landed"};

            std::string _current_state;

            // vector of all ForceTorques on vehicle
            std::vector<ForceTorque> _force_manager;

            std::vector<double> time;
            double mass;
            Vector3d cg_loc; // in body frame 

        public:

            void update();


    };
}