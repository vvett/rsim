#include "MassProperties.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

namespace py = pybind11;

PYBIND11_MODULE(_rsim, module) {
    module.doc() = "Python bindings for the rsim library";

    py::class_<rsim::Body>(module, "Body")
        .def(py::init<double, rsim::Body::Matrix3>(),
             py::arg("mass"), py::arg("moment_of_inertia"))
        .def_property("mass", &rsim::Body::mass, &rsim::Body::setMass,
                      "Body mass, which must be greater than zero.")
        .def_property(
            "moment_of_inertia",
            [](const rsim::Body& body) { return body.momentOfInertia(); },
            &rsim::Body::setMomentOfInertia,
            "The 3-by-3 moment-of-inertia matrix.")
        .def("print_mass", &rsim::Body::printMass)
        .def("__repr__", [](const rsim::Body& body) {
            return "rsim.Body(mass=" + std::to_string(body.mass()) + ")";
        });
}
