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
    Foam::functionObjects::theta

Group
    grpFieldFunctionObjects

Description
    Calculates the potential temperature, theta, from pressure and temperature

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    theta.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_theta_H
#define functionObjects_theta_H

#include "fieldExpression.H"
#include "dimensionedTypes.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class theta Declaration
\*---------------------------------------------------------------------------*/

class theta
:
    public fieldExpression
{
    // Private data
    
        //- Name of the pressure field, default is "p"
        word pName_;
        
        //- The reference pressure
        dimensionedScalar pRef_;
        
        //- R/Cp
        scalar kappa_;

    // Private Member Functions

        //- Calculate the vorticity field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("theta");


    // Constructors

        //- Construct from Time and dictionary
        theta
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~theta();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
