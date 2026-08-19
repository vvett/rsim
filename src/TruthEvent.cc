#include "TruthEvent.h"
#include "Parachute.h"
#include "Actuators.h"

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {
namespace {

enum class TokenKind { Identifier, Number, LParen, RParen, End };
struct Token {
    TokenKind kind;
    std::string text;
    double number = 0.0;
};

class Lexer {
public:
    explicit Lexer(const std::string& source) : source_(source) {}
    Token next() {
        while (position_ < source_.size() &&
               std::isspace(static_cast<unsigned char>(source_[position_]))) {
            ++position_;
        }
        if (position_ == source_.size()) {
            return {TokenKind::End, {}};
        }
        if (source_[position_] == '(') {
            ++position_;
            return {TokenKind::LParen, "("};
        }
        if (source_[position_] == ')') {
            ++position_;
            return {TokenKind::RParen, ")"};
        }
        const std::size_t begin = position_;
        if (std::isalpha(static_cast<unsigned char>(source_[position_])) ||
            source_[position_] == '_') {
            while (position_ < source_.size()) {
                const char c = source_[position_];
                if (!std::isalnum(static_cast<unsigned char>(c)) &&
                    c != '_' && c != '.') {
                    break;
                }
                ++position_;
            }
            return {TokenKind::Identifier, source_.substr(begin, position_ - begin)};
        }
        char* end = nullptr;
        const double number = std::strtod(source_.c_str() + position_, &end);
        if (end != source_.c_str() + position_) {
            position_ = static_cast<std::size_t>(end - source_.c_str());
            return {TokenKind::Number,
                    source_.substr(begin, position_ - begin), number};
        }
        throw std::invalid_argument("invalid character in truth-event expression");
    }
private:
    const std::string& source_;
    std::size_t position_ = 0;
};

double numeric(const StateValue& value) {
    return std::visit([](auto item) { return static_cast<double>(item); }, value);
}

}  // namespace

struct TruthEvent::Expression {
    enum class Kind { Comparison, And, Or } kind = Kind::Comparison;
    std::function<StateValue()> read;
    std::string operation;
    StateValue reference = 0.0;
    std::unique_ptr<Expression> left;
    std::unique_ptr<Expression> right;
    bool evaluate() const {
        if (kind == Kind::And) {
            return left->evaluate() && right->evaluate();
        }
        if (kind == Kind::Or) {
            return left->evaluate() || right->evaluate();
        }
        const double value = numeric(read());
        const double comparison = numeric(reference);
        if (!std::isfinite(value)) {
            return false;
        }
        if (operation == "greater_than") return value > comparison;
        if (operation == "less_than") return value < comparison;
        if (operation == "greater_than_or_equal") {
            return value >= comparison;
        }
        if (operation == "less_than_or_equal") {
            return value <= comparison;
        }
        if (operation == "equal_to") return value == comparison;
        return value != comparison;
    }
};

namespace {
class Parser {
public:
    Parser(const std::string& source,
           const std::unordered_map<std::string, std::function<StateValue()>>& values)
        : lexer_(source), values_(values), current_(lexer_.next()) {}
    std::unique_ptr<TruthEvent::Expression> parse() {
        auto result = parseOr();
        if (current_.kind != TokenKind::End) {
            throw std::invalid_argument(
                "unexpected token: " + current_.text);
        }
        return result;
    }
private:
    std::unique_ptr<TruthEvent::Expression> parseOr() {
        auto left = parseAnd();
        while (is("or")) {
            advance();
            auto node = std::make_unique<TruthEvent::Expression>();
            node->kind = TruthEvent::Expression::Kind::Or;
            node->left = std::move(left);
            node->right = parseAnd();
            left = std::move(node);
        }
        return left;
    }
    std::unique_ptr<TruthEvent::Expression> parseAnd() {
        auto left = parsePrimary();
        while (is("and")) {
            advance();
            auto node = std::make_unique<TruthEvent::Expression>();
            node->kind = TruthEvent::Expression::Kind::And;
            node->left = std::move(left);
            node->right = parsePrimary();
            left = std::move(node);
        }
        return left;
    }
    std::unique_ptr<TruthEvent::Expression> parsePrimary() {
        if (current_.kind == TokenKind::LParen) {
            advance();
            auto value = parseOr();
            require(TokenKind::RParen, ")");
            advance();
            return value;
        }
        require(TokenKind::Identifier, "state-variable alias");
        const std::string alias = current_.text;
        advance();
        require(TokenKind::Identifier, "comparison operator");
        const std::string operation = current_.text;
        advance();
        const bool recognized =
            operation == "greater_than" || operation == "less_than" ||
            operation == "greater_than_or_equal" ||
            operation == "less_than_or_equal" ||
            operation == "equal_to" || operation == "not_equal_to";
        if (!recognized) {
            throw std::invalid_argument(
                "unknown comparison operator: " + operation);
        }
        const auto found = values_.find(alias);
        if (found == values_.end()) {
            throw std::invalid_argument(
                "unknown state-variable alias: " + alias);
        }
        const StateValue current_value = found->second();
        const bool ordering = operation == "greater_than" ||
            operation == "less_than" ||
            operation == "greater_than_or_equal" ||
            operation == "less_than_or_equal";
        if (ordering && std::holds_alternative<bool>(current_value)) {
            throw std::invalid_argument("Boolean state variables support equality only: " + alias);
        }

        StateValue reference;
        if (current_.kind == TokenKind::Number) {
            if (std::holds_alternative<std::int64_t>(current_value)) {
                reference = static_cast<std::int64_t>(current_.number);
            } else if (std::holds_alternative<bool>(current_value)) {
                if (current_.number != 0.0 && current_.number != 1.0) {
                    throw std::invalid_argument(
                        "Boolean comparison literal must be true, false, 0, or 1");
                }
                reference = current_.number != 0.0;
            } else {
                reference = current_.number;
            }
        } else if (current_.kind == TokenKind::Identifier) {
            if (std::holds_alternative<bool>(current_value) &&
                (current_.text == "true" || current_.text == "false")) {
                reference = current_.text == "true";
            } else if (std::holds_alternative<std::int64_t>(current_value)) {
                static const std::unordered_map<std::string, std::int64_t> states{
                    {"ON_RAIL", 0}, {"THRUSTING", 1}, {"COASTING", 2},
                    {"DESCENDING", 3}, {"DROGUE_DEPLOYED", 4},
                    {"MAIN_DEPLOYED", 5}, {"LANDED", 6}};
                const auto state = states.find(current_.text);
                if (state == states.end()) {
                    throw std::invalid_argument(
                        "unknown symbolic state literal: " + current_.text);
                }
                reference = state->second;
            } else {
                throw std::invalid_argument(
                    "symbolic literal does not match state-variable type: " +
                    current_.text);
            }
        } else {
            throw std::invalid_argument("expected comparison literal");
        }
        auto result = std::make_unique<TruthEvent::Expression>();
        result->read = found->second;
        result->operation = operation;
        result->reference = reference;
        advance();
        return result;
    }
    bool is(const char* value) const {
        return current_.kind == TokenKind::Identifier && current_.text == value;
    }
    void require(TokenKind kind, const char* expected) {
        if (current_.kind != kind) {
            throw std::invalid_argument(std::string("expected ") + expected);
        }
    }
    void advance() { current_ = lexer_.next(); }
    Lexer lexer_;
    const std::unordered_map<std::string, std::function<StateValue()>>& values_;
    Token current_;
};
}  // namespace

TruthEvent::TruthEvent(std::string name, std::string trigger,
                       bool stop_simulation, bool one_shot)
    : name_(std::move(name)),
      trigger_(std::move(trigger)),
      stop_simulation_(stop_simulation),
      one_shot_(one_shot) {
    if (name_.empty() || trigger_.empty()) {
        throw std::invalid_argument(
            "truth event name and trigger must not be empty");
    }
}
TruthEvent::~TruthEvent() = default;
TruthEvent::TruthEvent(TruthEvent&&) noexcept = default;
TruthEvent& TruthEvent::operator=(TruthEvent&&) noexcept = default;
void TruthEvent::bind(
    const std::unordered_map<std::string,
                             std::function<StateValue()>>& values) {
    expression_ = Parser(trigger_, values).parse();
}
bool TruthEvent::evaluate() {
    if (one_shot_ && fired_) {
        return false;
    }
    if (!expression_ || !expression_->evaluate()) {
        return false;
    }
    fired_ = true;
    ++fire_count_;
    for (const auto& command : commands_) {
        command();
    }
    return true;
}
const std::string& TruthEvent::name() const noexcept { return name_; }
const std::string& TruthEvent::trigger() const noexcept { return trigger_; }
bool TruthEvent::fired() const noexcept { return fired_; }
std::uint64_t TruthEvent::fireCount() const noexcept { return fire_count_; }
void TruthEvent::addModelEnableCommand(Model& model, bool enabled) {
    commands_.push_back([&model, enabled]() { model.setEnabled(enabled); });
}
void TruthEvent::addParachuteDeployCommand(Parachute& parachute) {
    commands_.push_back([&parachute]() { parachute.deploy(); });
}
void TruthEvent::addActuatorCommand(
    IdealActuator& actuator, double value) {
    commands_.push_back(
        [&actuator, value]() { actuator.command(value); });
}

}  // namespace rsim
