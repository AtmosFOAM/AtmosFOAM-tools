/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InNamespace
    Foam::fvc

Description
    The divergence of a scalar over one time step

SourceFiles
    advect.C

\*---------------------------------------------------------------------------*/

#ifndef advect_H
#define advect_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "advectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    //- Calculate the divergence using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, fvPatchField, volMesh>> advect
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    );

    //- Calculate the divergence using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvPatchField, volMesh>> advect
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Calculate the divergence using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvPatchField, volMesh>> advect
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Calculate the divergence using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvPatchField, volMesh>> advect
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );

    //- Calculate the divergence using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvPatchField, volMesh>> advect
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "advect.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
