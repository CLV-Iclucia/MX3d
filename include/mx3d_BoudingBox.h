#ifndef MX3D_BOUNDINGBOX_H
#define MX3D_BOUNDINGBOX_H

#include "mx3d_Vector.h"
namespace mx3d
{
    struct BoundingBox
    {
        Vec3 lower;
        Vec3 upper;
        BoundingBox(const Vec3& _lower, const Vec3& _upper) : lower(_lower), upper(_upper) {}
        Vec3 size() const { return upper - lower; }
    };
}


#endif
