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
    Foam::fv::splitConvectionScheme

Description
    Abstract base class for a split convection using Gauss's divergence, first
    applying upwind then a higher-order correction

SourceFiles
    splitConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef splitConvectionScheme_H
#define splitConvectionScheme_H

#include "convectionScheme.H"
#include "upwind.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class splitConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class splitConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private Data
    
        tmp<gaussConvectionScheme<Type>> upwindConvection_;

public:

    // Constructors

        //- Construct from mesh and flux
        splitConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            upwindConvection_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        faceFlux,
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, faceFlux)
                        )
                    )
                )
            )
        {}


        //- Disallow default bitwise copy construction
        splitConvectionScheme(const splitConvectionScheme&) = delete;


    // Member Functions

        //- Acess to upwind convection scheme
        const gaussConvectionScheme<Type>& upwindConvect() const
        {
            return upwindConvection_();
        }

        //- The flux corrected high order correction to go from vf.oldTime()
        //- to vf via vfT which was updated using lowFlux
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        fluxCorrection
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& lowFlux,
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const GeometricField<Type, fvPatchField, volMesh>& vfT
        ) const = 0;
        
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const splitConvectionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "splitConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
