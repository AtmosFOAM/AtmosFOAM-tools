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
    Foam::fv::fluxCorrectedTransport

Description
    limits specified scheme with Zalesak FCT method, either to be monotonic or
    between specified limits.

SourceFiles
    fluxCorrectedTransport.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrectedTransport_H
#define fluxCorrectedTransport_H

#include "splitConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class fluxCorrectedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxCorrectedTransport
:
    public fv::splitConvectionScheme<Type>
{
    // Private Data

        //- Optional fixed limits (if both are zero then monotonic limiter)
        scalar Tmin_, Tmax_;

        //- convection scheme for the high order flux
        tmp<convectionScheme<Type>> tinterpHighScheme_;



public:

    //- Runtime type information
    TypeName("fluxCorrectedTransport");


    // Constructors

        //- Construct from components
        fluxCorrectedTransport
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            scalar Tmin,
            scalar Tmax,
            const tmp<convectionScheme<Type>>& tinterpHighScheme
        )
        :
            splitConvectionScheme<Type>(mesh, faceFlux),
            Tmin_(Tmin),
            Tmax_(Tmax),
            tinterpHighScheme_(tinterpHighScheme)
        {}

        //- Construct from mesh, faceFlux and Istream
        fluxCorrectedTransport
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            splitConvectionScheme<Type>(mesh, faceFlux),
            Tmin_(readScalar(is)),
            Tmax_(readScalar(is)),
            tinterpHighScheme_(convectionScheme<Type>::New(mesh, faceFlux, is))
        {}

        //- Disallow default bitwise copy construction
        fluxCorrectedTransport(const fluxCorrectedTransport&) = delete;


    // Member Functions
    
        //- The flux corrected high order correction to go from vf.oldTime()
        //- to vf via vfT which was updated using lowFlux
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxCorrection
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& lowFlux,
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const GeometricField<Type, fvPatchField, volMesh>& vfT
        ) const;
        
        const convectionScheme<Type>& interpHighScheme() const
        {
            return tinterpHighScheme_();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluxCorrectedTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluxCorrectedTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
