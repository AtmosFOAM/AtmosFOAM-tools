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

#include "pFromThetaExner.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pFromThetaExner, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        pFromThetaExner,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::pFromThetaExner::calc()
{
    if
    (
        
        foundObject<volScalarField>("Exner")
    )
    {
        Info<< "Reading thermophysical properties\n" << endl;
        
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                obr_.time().constant(),
                obr_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        const dimensionedScalar pRef
        (
            "pRef", dimPressure, readScalar(thermoDict.lookup("pRef"))
        );
        const dimensionedScalar T0
        (
            "T0", dimTemperature, readScalar(thermoDict.lookup("T0"))
        );
        const constTransport<hConstThermo<perfectGas<specie> > > air
        (
            thermoDict.subDict("mixture")
        );
        const dimensionedScalar R("R", dimGasConstant, air.R());
        const dimensionedScalar Cp
        (
            "Cp", dimGasConstant, air.Cp(pRef.value(),T0.value())
        );
        
        const volScalarField& Exner = lookupObject<volScalarField>("Exner");

        return store
        (
            resultName_,
            pRef*pow(Exner, Cp/R)
        );
    }
    else
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pFromThetaExner::pFromThetaExner
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "p")
{
    resultName_ = "p";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pFromThetaExner::~pFromThetaExner()
{}


// ************************************************************************* //
