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

#include "theta.H"
#include "addToRunTimeSelectionTable.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "constTransport.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(theta, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        theta,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::theta::calc()
{
    if
    (
        foundObject<volScalarField>(fieldName_)
     && foundObject<volScalarField>(pName_)
    )
    {
        const volScalarField& T = lookupObject<volScalarField>(fieldName_);
        const volScalarField& p = lookupObject<volScalarField>(pName_);

        return store
        (
            resultName_,
            T*pow(pRef_/p, kappa_)
        );
    }
    else
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::theta::theta
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "T"),
    pName_("p"),
    pRef_(dimensionedScalar("pRef", dimPressure, scalar(0))),
    kappa_(0)
{
    resultName_ = "theta";
    
    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    pRef_ = dimensionedScalar
    (
        "pRef", dimPressure, readScalar(thermoDict.lookup("pRef"))
    );
    dimensionedScalar T0
    (
        "T0", dimTemperature, readScalar(thermoDict.lookup("T0"))
    );
    const constTransport<hConstThermo<perfectGas<specie> > > air
    (
        thermoDict.subDict("mixture")
    );

    kappa_ = air.R()/air.Cp(pRef_.value(), T0.value());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::theta::~theta()
{}


// ************************************************************************* //
