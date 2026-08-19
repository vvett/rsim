#include "Kinematics.h"
#include "Simulator.h"

#include <Eigen/Dense>

#include <cassert>
#include <cmath>
#include <memory>

int main() {
    auto parent = std::make_shared<rsim::Frame>(
        "parent", rsim::InertialStatus::Inertial);
    rsim::ConstantForceTorque force_torque;
    rsim::ForceTorqueRegistry force_torques;
    force_torques.add(force_torque);
    rsim::RigidBody rigid_body(10.0, rsim::Matrix3d::Identity());

    rsim::Kinematics first("body-1", parent, force_torques, rigid_body);
    rsim::Kinematics second("body-2", parent, force_torques, rigid_body);

    rsim::FrameGraph frames;
    frames.addFrame(parent);
    rsim::Kinematics named(
        "body-named", frames, "parent", force_torques, rigid_body);

    assert(first.bodyFrame);
    assert(first.bodyFrameDefinition);
    assert(first.bodyFrameMotion);
    assert(first.bodyFrame->parent() == parent);
    assert(!first.bodyFrame->isInertial());
    assert(first.bodyFrameMotion ==
           &first.bodyFrameDefinition->motionIntoParent());
    assert(first.bodyFrameMotion == &first.bodyFrame->motionIntoParent());
    assert(first.bodyFrameMotion != second.bodyFrameMotion);
    assert(named.integrationFrame() == parent);

    first.update();
    assert(first.bodyFrameMotion->transform().translation().isZero());
    assert(first.bodyFrameMotion->transform().rotation().isApprox(
        Eigen::Quaterniond::Identity()));
    assert(first.bodyFrameMotion->childOriginVelocity().isZero());
    assert(first.bodyFrameMotion->childOriginAcceleration().isZero());
    assert(first.bodyFrameMotion->childAngularVelocity().isZero());
    assert(first.bodyFrameMotion->childAngularAcceleration().isZero());

    // A constant body-x force advances position and velocity in the same RK4 step.
    force_torque.setLoad(rsim::Vector3d{20.0, 0.0, 0.0},
                         rsim::Vector3d::Zero());
    rsim::Simulator simulator(0.1);
    simulator.addModel(first);
    simulator.step();
    assert(first.velocityInParent().isApprox(
        rsim::Vector3d{0.2, 0.0, 0.0}, 1e-12));
    assert(first.positionInParent().isApprox(
        rsim::Vector3d{0.01, 0.0, 0.0}, 1e-12));

    // Constant angular velocity rotates attitude during the coupled RK4 step.
    force_torque.setLoad(rsim::Vector3d::Zero(), rsim::Vector3d::Zero());
    first.setState(rsim::Vector3d::Zero(), rsim::Vector3d::Zero(),
                   Eigen::Quaterniond::Identity(),
                   rsim::Vector3d{0.0, 0.0, 1.0});
    simulator.step();
    const Eigen::Quaterniond expected_rotation{
        Eigen::AngleAxisd(0.1, rsim::Vector3d::UnitZ())};
    assert(first.rotationIntoParent().isApprox(expected_rotation, 1e-7));

    // A force through an offset CG translates the body without rotating it.
    rsim::ConstantForceTorque cg_force("force-through-cg");
    rsim::ForceTorqueRegistry cg_forces("offset-cg-forces");
    cg_forces.add(cg_force);
    const rsim::Vector3d offset_cg{2.0, 0.0, 0.0};
    rsim::RigidBody offset_body(
        10.0, rsim::Matrix3d::Identity(), offset_cg, "offset-body");
    rsim::Kinematics offset_kinematics(
        "offset-kinematics", parent, cg_forces, offset_body);
    cg_force.setLoad(
        rsim::Vector3d{0.0, 10.0, 0.0}, rsim::Vector3d::Zero(), offset_cg);
    rsim::Simulator offset_simulator(0.1);
    offset_simulator.addModel(offset_kinematics);
    offset_simulator.step();
    assert(offset_kinematics.angularVelocityInBody().isZero(1e-12));
    assert(offset_kinematics.velocityInParent().isApprox(
        rsim::Vector3d{0.0, 0.1, 0.0}, 1e-12));
}
