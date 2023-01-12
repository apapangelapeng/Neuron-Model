#if !defined(STRING_PERCISION_H)
#define STRING_PERCISION_H 1


#include <string.h>
#include <sstream>
template <typename T>


std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


#endif
