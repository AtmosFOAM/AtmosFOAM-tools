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

Class
    Foam::functionObjects::turbulentLength

Group
    grpFieldFunctionObjects

Description
    Evaluates and writes the turbulence length scale, k^1.5/epsilon

SourceFiles
    turbulentLength.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulentLength_H
#define functionObjects_turbulentLength_H

#include "fieldExpression.H"
#include "dimensionedTypes.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class turbulentLength Declaration
\*---------------------------------------------------------------------------*/

class turbulentLength
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the length scale and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("turbulentLength");


    // Constructors

        //- Construct from Time and dictionary
        turbulentLength
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulentLength();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
