#include "phist_version.hpp"

namespace phist
{
    extern "C"
    {
        std::string version()
        {
            return std::string{"phist version "};
        }
    }

}