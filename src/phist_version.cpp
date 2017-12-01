#include "phist_version.hpp"

namespace phist
{
    extern "C"
    {
        constexpr std::string version()
        {
            return std::string{"phist version "};
        }
    }

}