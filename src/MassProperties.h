#pragma once
#include <array>

namespace rsim {

class Body {
public:
    using Matrix3 = std::array<std::array<double, 3>, 3>;

    Body(double mass, Matrix3 moment_of_inertia);

    [[nodiscard]] double mass() const noexcept;
    void setMass(double mass);

    [[nodiscard]] const Matrix3& momentOfInertia() const noexcept;
    void setMomentOfInertia(Matrix3 moment_of_inertia) noexcept;

    void printMass() const;

private:
    Matrix3 moment_of_inertia_;
    double mass_;
};

}  // namespace rsim
