#ifndef MX3D_BOUNDINGBOX_H
#define MX3D_BOUNDINGBOX_H

#include "mx3d_Vector.h"
namespace mx3d
{
    template<uint Dim>
    struct BoundingBox
    {
        Vector<Real, Dim> lower;
        Vector<Real, Dim> upper;
        BoundingBox(Vector<Real, Dim> _lower, Vector<Real, Dim> _upper) : lower(std::move(_lower)), upper(std::move(_upper)) {}
        Vector<Real, Dim> size() const { return upper - lower; }
    };
}


#endif
