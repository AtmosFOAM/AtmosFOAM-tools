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
    Foam::fv::MPDATA_theta

Description
    MPDATA_theta is MPDATA with Crank-Nicholson time stepping.

SourceFiles
    MPDATA_theta.C

\*---------------------------------------------------------------------------*/

#ifndef MPDATA_theta_H
#define MPDATA_theta_H

#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class MPDATA_theta Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MPDATA_theta
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Dictionary for reading the MPDATA parameters
        dictionary dict_;

        //- Number of corrections (1 is normal)
        label nCorr_;
        
        // MPDATA gauge (the minimum allowed value)
        const scalar gauge_;

        //- ante-diffusive flux
        mutable surfaceScalarField* anteDPtr_;

        //- upwind interpolation
        tmp<gaussConvectionScheme<Type>> tupwindConvection_;
        
        //- upwind interpolation using ante-diffusive flux
        tmp<gaussConvectionScheme<Type>> tanteDConvect_;

public:

    //- Runtime type information
    TypeName("MPDATA_theta");

    // Constructors

        //- Construct from components
        MPDATA_theta
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            label nCorr__ = 1,
            scalar gauge__ = 0
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            dict_(),
            nCorr_(nCorr__),
            gauge_(gauge__),
            anteDPtr_(new surfaceScalarField(faceFlux)),
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
            tanteDConvect_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        anteD(),
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, anteD())
                        )
                    )
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        MPDATA_theta
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            dict_(is),
            nCorr_(readLabel(dict_.lookup("nCorr"))),
            gauge_(readScalar(dict_.lookup("gauge"))),
            anteDPtr_(new surfaceScalarField(faceFlux)),
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
            tanteDConvect_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        anteD(),
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, anteD())
                        )
                    )
                )
            )
        {}
        
        //- Disallow default bitwise copy construction
        MPDATA_theta(const MPDATA_theta&) = delete;

        //- Destructor
        ~MPDATA_theta()
        {
            delete anteDPtr_;
        }

    // Member Functions
    
        // Return the ante-diffusive flux
        surfaceScalarField& anteD() const { return *anteDPtr_; }
    
        //- Access to upwind convection scheme
        const gaussConvectionScheme<Type>& upwindConvect() const
        {
            return tupwindConvection_();
        }

        //- Access to anti-diffusive convection scheme
        const gaussConvectionScheme<Type>& anteDConvect() const
        {
            return tanteDConvect_();
        }
        
        //- Calculate the anti-diffusive flux
        void calculateAnteD
        (
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& offCentre
        ) const;

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

        // The part of the advection which is solved implicitly
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
        void operator=(const MPDATA_theta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MPDATA_theta.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
