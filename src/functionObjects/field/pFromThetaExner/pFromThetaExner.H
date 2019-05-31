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

Class
    Foam::functionObjects::pFromThetaExner

Group
    grpFieldFunctionObjects

Description
    Calculates the pressure p from Exner, pRef, cP, R via the relation:
        p = pRef*(Exner)^(cP/R).

SourceFiles
    pFromThetaExner.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pFromThetaExner_H
#define functionObjects_pFromThetaExner_H

#include "fieldExpression.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "OFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "constTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class pFromThetaExner Declaration
\*---------------------------------------------------------------------------*/

class pFromThetaExner
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate p from Exner
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pFromThetaExner");


    // Constructors

        //- Construct from Time and dictionary
        pFromThetaExner
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pFromThetaExner();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
