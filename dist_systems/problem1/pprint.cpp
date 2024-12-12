// Based on
// https://stackoverflow.com/questions/70905100/c-concepts-for-numeric-types-and-preventing-instantiation-on-pointers
// https://www.geeksforgeeks.org/how-to-iterate-over-the-elements-of-an-stdtuple-in-c/
#include <iostream>
#include <sstream>
#include <type_traits>
#include <string>
#include <set>
#include <vector>
#include <utility>
#include <tuple>


template <typename T>
concept IntegerType = requires(T param) {
    requires std::is_integral_v<T>;
};

struct TPrettyPrinter {

    TPrettyPrinter() : str{} {}

    std::string Str() const {
        return str;
    }

    template <IntegerType T>
    TPrettyPrinter& f(const T& value) {
        str += std::to_string(value);
        return *this;
    }

    TPrettyPrinter& f(const std::string& value) {
        str += value;
        return *this;
    }

    template<typename T>
    TPrettyPrinter& f(const std::vector<T>& value) {
        return format_container(value, "[", "]");
    }

    template<typename T>
    TPrettyPrinter& f(const std::set<T>& value) {
        return format_container(value, "{", "}");
    }

    template<typename... Args>
    TPrettyPrinter& f(const std::tuple<Args...>& value) {
        size_t length = sizeof...(Args);
        std::stringstream ss{};
        std::apply(
            [length, &ss](auto const&... ps) {
                ss << "(";
                size_t k = 0;
                ((ss << TPrettyPrinter().f(ps).Str() << (++k == length ? "" : ", ")), ...);
                ss << ")";
            },
            value);
        str += ss.str();
        return *this;
    }

    template<typename T1, typename T2>
    TPrettyPrinter& f(const std::pair<T1, T2>& value) {
        str += "<" + TPrettyPrinter().f(value.first).f(", ").f(value.second).Str() + ">";
        return *this;
    }

private:
    std::string str;

    template<typename C>
    TPrettyPrinter& format_container(
            const C& value,
            std::string open_bracket,
            std::string close_bracket) {
        // str += open_bracket + TPrettyPrinter().f(*(value.begin())).Str();
        str += open_bracket;
        for (auto it = value.begin(); it != value.end(); it++) {
            str += TPrettyPrinter().f(*it).Str() + ", ";
        }
        str = str.substr(0, str.size() - 2);
        str += close_bracket;
        return *this;
    }
};


template<typename T>
void pprint(const T& t) {
    std::cout << TPrettyPrinter().f(t).Str() << std::endl;
}

template<typename T>
std::string spprint(const T& t) {
    return TPrettyPrinter().f(t).Str();
}
