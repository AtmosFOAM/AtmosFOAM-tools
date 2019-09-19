/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "convection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(convection, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        convection,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::convection::calc()
{
    if
    (
        foundObject<volScalarField>(fieldNames_[0])
     && foundObject<volVectorField>(fieldNames_[1])
    )
    {
        const volScalarField& theta(lookupObject<volScalarField>(fieldNames_[0]));
        const volVectorField& u(lookupObject<volVectorField>(fieldNames_[1]));

        volVectorField gradTheta = fvc::grad(theta);
        volScalarField dThetadz = gradTheta.component(2);
        volScalarField w = u.component(2);

        tmp<volScalarField> wdThetadz(-w*dThetadz);
        wdThetadz.ref() *= 0.25*(sign(w)+1) * (1-sign(dThetadz));
        return store(resultName_, wdThetadz);
    }
    else
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::convection::convection
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldsExpression(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::convection::~convection()
{}


// ************************************************************************* //
