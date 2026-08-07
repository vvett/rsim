#pragma once
#include <array>
#include <Eigen/Dense>
#include "defs.h"

namespace rsim {
    
    class RigidBody {

        private:
            Matrix3d moment_of_inertia_;
            double mass_;

        public:
            RigidBody(double mass, Matrix3d moment_of_inertia);

            double mass() const;

            void setMass(double mass);

            const Matrix3d& momentOfInertia() const noexcept;
            
            void setMomentOfInertia(Matrix3d moment_of_inertia);

            void printMass() const;

    };
}  // namespace rsim
