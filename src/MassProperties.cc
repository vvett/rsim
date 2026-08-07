#include "MassProperties.h"
#include <iostream>
#include <stdexcept>

namespace rsim {

RigidBody::RigidBody(double mass, Matrix3d moment_of_inertia)
    : moment_of_inertia_(moment_of_inertia), mass_(mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
}

double RigidBody::mass() const {
    return mass_;
}

void RigidBody::setMass(double mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
    mass_ = mass;
}

const Matrix3d& RigidBody::momentOfInertia() const noexcept {
    return moment_of_inertia_;
}

void RigidBody::setMomentOfInertia(Matrix3d moment_of_inertia) {
    moment_of_inertia_ = moment_of_inertia;
}

void RigidBody::printMass() const {
    std::cout << "Body mass is " << mass_ << '\n';
}

}  // namespace rsim
