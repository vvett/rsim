#include "MassProperties.h"
#include <iostream>
#include <stdexcept>

namespace rsim {

Body::Body(double mass, Matrix3 moment_of_inertia)
    : moment_of_inertia_(moment_of_inertia), mass_(mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
}

double Body::mass() const noexcept {
    return mass_;
}

void Body::setMass(double mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
    mass_ = mass;
}

const Body::Matrix3& Body::momentOfInertia() const noexcept {
    return moment_of_inertia_;
}

void Body::setMomentOfInertia(Matrix3 moment_of_inertia) noexcept {
    moment_of_inertia_ = moment_of_inertia;
}

void Body::printMass() const {
    std::cout << "Body mass is " << mass_ << '\n';
}

}  // namespace rsim
