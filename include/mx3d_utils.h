//
// Created by creeper on 23-3-28.
//

#ifndef MX3D_MX3D_UTILS_H
#define MX3D_MX3D_UTILS_H
#include "mx3d_types.h"
namespace mx3d
{
    static const Real EPS = 1e-14;
    template<typename T>
    inline bool isZero(T x)
    {
        if constexpr (is_same_v<T, float> || is_same_v<T, Real>)
            return x <= EPS && x >= -EPS;;
        else return x == 0;
    }
    template<typename T>
    inline bool isEqual(T x, T y)
    {
        if constexpr (is_same_v<T, float> || is_same_v<T, Real>)
        return x - y <= EPS && x - y >= -EPS;;
        else return x == y;
    }
}
#endif //MX3D_MX3D_UTILS_H
