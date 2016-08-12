#ifndef sphericalVelocity_H
#define sphericalVelocity_H

#include "dimensionedScalar.H"
#include "dimensionedVector.H"

namespace Foam
{

class sphericalVelocity
{
    public:
    sphericalVelocity(dimensionedScalar u, dimensionedScalar v, dimensionedScalar w);
    dimensionedVector toCartesian(void) const;

    private:
    const dimensionedVector v;
};

}

#endif
