#include "utility/MatFileReader.h"

#include <matio.h>

#include <numeric>
#include <stdexcept>
#include <utility>

namespace rsim::utility {

namespace {

// RAII releases variables returned by matio even when validation throws.
struct MatVariableCloser {
    void operator()(matvar_t* variable) const noexcept {
        Mat_VarFree(variable);
    }
};

using MatVariable = std::unique_ptr<matvar_t, MatVariableCloser>;

std::size_t elementCount(const matvar_t& variable) {
    return std::accumulate(
        variable.dims,
        variable.dims + variable.rank,
        std::size_t{1},
        std::multiplies<>());
}

}  // namespace

class MatFileReader::Implementation {
public:
    explicit Implementation(const std::string& path)
        : file_(Mat_Open(path.c_str(), MAT_ACC_RDONLY)) {
        if (file_ == nullptr) {
            throw std::runtime_error("cannot open MAT file: " + path);
        }
    }

    ~Implementation() { Mat_Close(file_); }

    // Read one required variable and transfer its allocation to a smart pointer.
    [[nodiscard]] MatVariable read(const std::string& name) const {
        MatVariable variable(Mat_VarRead(file_, name.c_str()));
        if (!variable) {
            throw std::runtime_error("MAT file is missing " + name);
        }
        if (variable->isComplex || variable->data == nullptr) {
            throw std::runtime_error("MAT variable must be real: " + name);
        }
        return variable;
    }

private:
    mat_t* file_;
};

MatFileReader::MatFileReader(const std::string& path)
    : implementation_(std::make_unique<Implementation>(path)) {}

MatFileReader::~MatFileReader() = default;
MatFileReader::MatFileReader(MatFileReader&&) noexcept = default;
MatFileReader& MatFileReader::operator=(MatFileReader&&) noexcept = default;

MatDoubleArray MatFileReader::readDoubleArray(const std::string& name) const {
    const MatVariable variable = implementation_->read(name);
    if (variable->class_type != MAT_C_DOUBLE ||
        variable->data_type != MAT_T_DOUBLE) {
        throw std::runtime_error("MAT variable must contain doubles: " + name);
    }
    const auto* begin = static_cast<const double*>(variable->data);
    return {
        {variable->dims, variable->dims + variable->rank},
        {begin, begin + elementCount(*variable)}
    };
}

std::uint32_t MatFileReader::readUint32(const std::string& name) const {
    const MatVariable variable = implementation_->read(name);
    if (variable->class_type != MAT_C_UINT32 ||
        variable->data_type != MAT_T_UINT32 ||
        elementCount(*variable) != 1U) {
        throw std::runtime_error("MAT variable must be a uint32 scalar: " + name);
    }
    return *static_cast<const std::uint32_t*>(variable->data);
}

std::string MatFileReader::readUtf8(const std::string& name) const {
    const MatVariable variable = implementation_->read(name);
    if (variable->class_type != MAT_C_UINT8 ||
        variable->data_type != MAT_T_UINT8) {
        throw std::runtime_error("MAT variable must contain UTF-8 bytes: " + name);
    }
    const auto* begin = static_cast<const char*>(variable->data);
    return {begin, begin + elementCount(*variable)};
}

}  // namespace rsim::utility
