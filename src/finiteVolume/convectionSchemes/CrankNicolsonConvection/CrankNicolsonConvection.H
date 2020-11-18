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
    Foam::fv::CrankNicolsonConvection

Description
    Applies linear differencing with Crank-Nicolson time stepping as an explicit
    correction on upwind and Backward Euler

SourceFiles
    CrankNicolsonConvection.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolsonConvection_H
#define CrankNicolsonConvection_H

#include "upwind.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class CrankNicolsonConvection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CrankNicolsonConvection
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Off Centering
        scalar offCenter_;
        
        //- Number of high order corrections
        label nCorr_;
        
        //- upwind interpolation
        tmp<gaussConvectionScheme<Type>> tupwindConvection_;

        //- convection scheme for the high order flux
        tmp<gaussConvectionScheme<Type>> thighConvection_;

public:

    //- Runtime type information
    TypeName("CrankNicolsonConvection");


    // Constructors

        //- Construct from components
        CrankNicolsonConvection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            scalar offCenter,
            label nCorr,
            const tmp<surfaceInterpolationScheme<Type>>& tinterpHighScheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            offCenter_(offCenter),
            nCorr_(nCorr),
            tupwindConvection_
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
            ),
            thighConvection_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        faceFlux,
                        tinterpHighScheme
                    )
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        CrankNicolsonConvection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            offCenter_(readScalar(is)),
            nCorr_(readLabel(is)),
            tupwindConvection_
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
            ),
            thighConvection_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        faceFlux,
                        surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
                    )
                )
            )
        {}

        //- Disallow default bitwise copy construction
        CrankNicolsonConvection(const CrankNicolsonConvection&) = delete;


    // Member Functions
    
        //- Off centering coefficient
        scalar offCenter() { return offCenter_; }
        
        //- Number of high order correction steps >=1
        label nCorr() { return nCorr_;}

        //- Access to upwind convection scheme
        const gaussConvectionScheme<Type>& upwindConvect() const
        {
            return tupwindConvection_();
        }

        const gaussConvectionScheme<Type>& highConvection() const
        {
            return thighConvection_();
        }

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
        void operator=(const CrankNicolsonConvection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CrankNicolsonConvection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
