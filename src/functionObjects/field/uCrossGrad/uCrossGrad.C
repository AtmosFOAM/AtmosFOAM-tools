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

#include "uCrossGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(uCrossGrad, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        uCrossGrad,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::uCrossGrad::calc()
{
    if
    (
        foundObject<volScalarField>(fieldNames_[0])
     && foundObject<volVectorField>(fieldNames_[1])
    )
    {
        const volScalarField& T(lookupObject<volScalarField>(fieldNames_[0]));
        const volVectorField& u
        (
            lookupObject<volVectorField>(fieldNames_[1])
        );

        tmp<volVectorField> uxgradT(u ^ fvc::grad(T));
        return store(resultName_, uxgradT);
    }
    else
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::uCrossGrad::uCrossGrad
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

Foam::functionObjects::uCrossGrad::~uCrossGrad()
{}


// ************************************************************************* //
