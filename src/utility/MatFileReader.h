#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace rsim::utility {

// Keeps a numeric array's MATLAB dimensions next to its column-major values.
struct MatDoubleArray {
    std::vector<std::size_t> dimensions;
    std::vector<double> values;
};

// Provides the small typed MAT-file interface needed by rsim databases.
class MatFileReader {
public:
    explicit MatFileReader(const std::string& path);
    ~MatFileReader();

    MatFileReader(const MatFileReader&) = delete;
    MatFileReader& operator=(const MatFileReader&) = delete;
    MatFileReader(MatFileReader&&) noexcept;
    MatFileReader& operator=(MatFileReader&&) noexcept;

    // Read one real double array and preserve its MATLAB dimensions.
    [[nodiscard]] MatDoubleArray readDoubleArray(const std::string& name) const;

    // Read one scalar stored with MATLAB's uint32 type.
    [[nodiscard]] std::uint32_t readUint32(const std::string& name) const;

    // Read a uint8 array as a UTF-8 byte string.
    [[nodiscard]] std::string readUtf8(const std::string& name) const;

private:
    // The implementation hides matio types from every public rsim header.
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace rsim::utility
