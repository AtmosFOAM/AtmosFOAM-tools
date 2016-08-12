#include "sphericalVelocity.H"

#include "vector.H"

Foam::sphericalVelocity::sphericalVelocity
(
        dimensionedScalar u, dimensionedScalar v, dimensionedScalar w
)
:
    v("sphericalVelocity", dimVelocity, vector(u.value(),v.value(),w.value()))
{};

Foam::dimensionedVector Foam::sphericalVelocity::toCartesian() const
{
    return v;
}

