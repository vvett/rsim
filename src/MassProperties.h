#pragma once
#include <array>
#include <Eigen/Dense>
#include "defs.h"
#include <vector>
#include "Frames.h"

namespace rsim {
    
    class RigidBody {

        private:

            // total moment of inertia
            Matrix3d moment_of_inertia_;

            // dry moment of inertia (cg)
            Matrix3d dry_moi;

            // total mass
            double mass_;

            // dry mass
            double dry_mass_;

            // Body frame
            std::vector<double> time;

            // Body frame
            Frame* body_frame;

        public:

            RigidBody(double mass, Matrix3d moment_of_inertia);

            double mass() const;

            void setMass(double mass);

            Matrix3d momentOfInertia() const noexcept;
            
            void setMomentOfInertia(Matrix3d moment_of_inertia);

            void printMass() const;

            Vector3d cg() const;

            void update();

    };
}  // namespace rsim
