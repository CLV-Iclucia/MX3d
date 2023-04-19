#ifndef MX3D_BOUNDINGBOX_H
#define MX3D_BOUNDINGBOX_H

#include "mx3d_Vector.h"
namespace mx3d
{
    template<uint Dim>
    struct BoundingBox
    {
        static_assert(Dim == 2u || Dim == 3u, "Only 2D and 3D bounding boxes are supported.");
        Vector<Real, Dim> lower;
        Vector<Real, Dim> upper;
        BoundingBox(Vector<Real, Dim> _lower, Vector<Real, Dim> _upper) : lower(std::move(_lower)), upper(std::move(_upper)) {}
        Vector<Real, Dim> size() const { return upper - lower; }
        BoundingBox merge(const BoundingBox& bbox) const
        {
            if constexpr (Dim == 2)
                return {{std::min(lower.x, bbox.lower.x), std::max(upper.x, bbox.upper.x)},
                        {std::min(lower.y, bbox.lower.y), std::max(upper.y, bbox.upper.y)}};
            else if constexpr (Dim == 3)
                return {{std::min(lower.x, bbox.lower.x), std::max(upper.x, bbox.upper.x)},
                        {std::min(lower.y, bbox.lower.y), std::max(upper.y, bbox.upper.y)},
                        {std::min(lower.z, bbox.lower.z), std::max(upper.z, bbox.upper.z)}};
        }
    };
}


#endif
