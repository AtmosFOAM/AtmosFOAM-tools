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
    Foam::functionObjects::uCrossGrad

Group
    grpFieldFunctionObjects

Description
    Calculates the u x grad(T)

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    uCrossGrad.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_uCrossGrad_H
#define functionObjects_uCrossGrad_H

#include "fieldsExpression.H"
#include "dimensionedTypes.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class uCrossGrad Declaration
\*---------------------------------------------------------------------------*/

class uCrossGrad
:
    public fieldsExpression
{
    // Private Member Functions

        //- Calculate the vorticity field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("uCrossGrad");


    // Constructors

        //- Construct from Time and dictionary
        uCrossGrad
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uCrossGrad();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
