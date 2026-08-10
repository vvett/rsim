#include "MassProperties.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <string>

namespace py = pybind11;

PYBIND11_MODULE(_rsim, module) {
    module.doc() = "Python bindings for the rsim library";

    py::class_<rsim::RigidBody>(module, "RigidBody")
        .def(py::init<double, rsim::Matrix3d>(),
             py::arg("mass"), py::arg("moment_of_inertia"))
        .def_property("mass", &rsim::RigidBody::mass, &rsim::RigidBody::setMass,
                      "Body mass, which must be greater than zero.")
        .def_property(
            "moment_of_inertia",
            [](const rsim::RigidBody& rigidBody) { return rigidBody.momentOfInertia(); },
            &rsim::RigidBody::setMomentOfInertia,
            "The 3-by-3 moment-of-inertia matrix.")
        .def("print_mass", &rsim::RigidBody::printMass)
        .def("__repr__", [](const rsim::RigidBody& rigidBody) {
            return "rsim.RigidBody(mass=" + std::to_string(rigidBody.mass()) + ")";
        });
}
