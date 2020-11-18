/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::advectionScheme

Description
    Abstract base class for advection schemes.

SourceFiles
    advectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef advectionScheme_H
#define advectionScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                 Class advectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class advectionScheme
:
    public tmp<advectionScheme<Type>>::refCount
{
    // Private Data

        //- Hold reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("advectionScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            advectionScheme,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh
        advectionScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        advectionScheme(const advectionScheme&);



    // Selectors

        //- Return new tmp advection scheme
        static tmp<advectionScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~advectionScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }


        //- Return the divergence of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, fvPatchField, volMesh>>
        advect
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const surfaceScalarField& faceFlux
        ) const;

        //- Return the divergence of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        advect
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            const surfaceScalarField& faceFlux
        ) const;

        //- Return the divergence of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, fvPatchField, volMesh>>
        advect
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<surfaceScalarField>& tFaceFlux
        ) const;

        //- Return the divergence of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        advect
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            const tmp<surfaceScalarField>& tFaceFlux
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const advectionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeadvectionTypeScheme(SS, Type)                                      \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
advectionScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>                 \
    add##SS##Type##MeshFluxConstructorToTable_;

#define makeadvectionScheme(SS)                                                \
                                                                               \
makeadvectionTypeScheme(SS, scalar)                                            \
makeadvectionTypeScheme(SS, vector)                                            \
makeadvectionTypeScheme(SS, sphericalTensor)                                   \
makeadvectionTypeScheme(SS, symmTensor)                                        \
makeadvectionTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "advectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
