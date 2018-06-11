/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "frictionVelocity.H"
#include "turbulenceModel.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(frictionVelocity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        frictionVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::frictionVelocity::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        tmp<volScalarField> nuEff
        (
            mesh_.lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            ).nuEff()
        );
        
        const volVectorField& U = mesh_.lookupObject<volVectorField>(fieldName_);

        return store
        (
            resultName_,
            sqrt(fvc::interpolate(nuEff)*mag(fvc::snGrad(mag(U))))
        );
    }
    else
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::frictionVelocity::frictionVelocity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName("frictionVelocity", "U");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::frictionVelocity::~frictionVelocity()
{}


// ************************************************************************* //
