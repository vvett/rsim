#pragma once

#include "Frames.h"

namespace rsim {

/** A Julian date split into two doubles to retain sub-second precision. */
struct JulianDate {
    double first_part = 0.0;
    double second_part = 0.0;
};

/**
 * A simple ECEF definition relative to fixed J2000 axes.
 *
 * The supplied date is treated directly as mean universal time. The model uses
 * Greenwich mean sidereal rotation and deliberately omits UT1 corrections,
 * precession, nutation, polar motion, and relativistic time scales.
 */
class J2000EarthRotation final : public FrameDefinition {
public:
    /**
     * @param current_date Live Julian date owned and advanced by the caller.
     *        It must outlive this definition.
     */
    explicit J2000EarthRotation(const JulianDate& current_date);

    void update() override;
    [[nodiscard]] const FrameMotion& motionIntoParent() const override;

private:
    const JulianDate& current_date_;
    FrameMotion j2000_from_ecef_;
};

}  // namespace rsim
