/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gaussian.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussian::gaussian()
:
    centre_(dimensioned<point>("centre", dimLength, point::zero)),
    radius_(dimensionedScalar("radius", dimLength, scalar(0))),
    max_(dimensionedScalar("T", dimless, scalar(0))),
    periodMin_(point::zero),
    periodMax_(point::zero)
{}


Foam::gaussian::gaussian
(
    const dimensioned<point>& centre__,
    const dimensionedScalar radius__,
    const dimensionedScalar max__,
    const point& periodMin__,
    const point& periodMax__
)
:
    centre_(centre__),
    radius_(radius__),
    max_(max__),
    periodMin_(periodMin__),
    periodMax_(periodMax__)
{}


Foam::gaussian::gaussian(const gaussian& g)
:
    centre_(g.centre()),
    radius_(g.radius()),
    max_(g.max()),
    periodMin_(g.periodMin()),
    periodMax_(g.periodMax())
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussian::~gaussian()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::gaussian::field(const fvMesh& mesh) const
{
    tmp<volScalarField> tvf
    (
        new volScalarField
        (
            IOobject(max_.name(), mesh.time().timeName(), mesh),
            mesh,
            dimensionedScalar(max_.name(), max_.dimensions(), scalar(0))
        )
    );
    volScalarField& vf = tvf.ref();
    
    volVectorField dist = mesh.C() - centre_;
    
    vector period = periodMax_ - periodMin_;
    for(label dim = 0; dim < 3; dim++)
    {
        if (period[dim] < VSMALL) period[dim] = VGREAT;
    }
    forAll(dist, cellI)
    {
        for(label dim = 0; dim < 3; dim++)
        {
            if (dist[cellI][dim] > 0.5*period[dim])
            {
                dist[cellI][dim] -= period[dim];
            }
            else if (dist[cellI][dim] < -0.5*period[dim])
            {
                dist[cellI][dim] += period[dim];
            }
        }
    }
    
    vf = max_*exp(-0.5*magSqr(dist)/sqr(radius_));
    
    return tvf;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::gaussian::operator=(const gaussian& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::gaussian::operator=(const Foam::gaussian&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
    centre_ = rhs.centre();
    radius_ = rhs.radius();
    max_ = rhs.max();
}


// ************************************************************************* //
