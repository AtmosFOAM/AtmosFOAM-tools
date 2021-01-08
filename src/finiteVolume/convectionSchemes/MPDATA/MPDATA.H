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
    Foam::fv::MPDATA

Description
    MPDATA advection either as an explicit or implicit correction on first order
    forward or backward time stepping. With or without a guage. 

SourceFiles
    MPDATA.C

\*---------------------------------------------------------------------------*/

#ifndef MPDATA_H
#define MPDATA_H

#include "upwind.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class MPDATA Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MPDATA
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Pointer to the MPDATA gauge (minimum allowed value)
        mutable dimensionedScalar* gaugePtr_;
        
        //- Maximum Courant number for explicit advection
        scalar maxCoExp_;
        
        //- ante-diffusive flux
        mutable surfaceScalarField* anteDPtr_;

        //- upwind interpolation
        tmp<gaussConvectionScheme<Type>> tupwindConvection_;
        
        //- upwind interpolation using ante-diffusive flux
        tmp<gaussConvectionScheme<Type>> tanteDConvect_;

public:

    //- Runtime type information
    TypeName("MPDATA");

    //- Flags for implicit or explicit
    enum implicitExplicit
    {
        IMPLICIT,
        EXPLICIT
    };

    // Constructors

        //- Construct from components
        MPDATA
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            scalar gauge__ = 0,
            scalar maxCoExp__ = 1
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            gaugePtr_(new dimensionedScalar("gauge", dimless, gauge__)),
            maxCoExp_(maxCoExp__),
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
        MPDATA
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            gaugePtr_(new dimensionedScalar("gauge", dimless, readScalar(is))),
            maxCoExp_(readScalar(is)),
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
        {
            if (gauge().value() < 0)
            {
                FatalErrorIn("MPDATA::MPDATA") << " gauge must be >=0 but "
                    << gauge().value() << " was given\n" << exit(FatalError);
            }
            if (maxCoExp_ < 0)
            {
                FatalErrorIn("MPDATA::MPDATA") << " maxCoExp must be >=0 but "
                    << maxCoExp_ << " was given\n" << exit(FatalError);
            }
        }

        //- Disallow default bitwise copy construction
        MPDATA(const MPDATA&) = delete;


    // Member Functions
    
        // The gauge as a dimensionedScalar (the MPDATA minimum allowed value)
        dimensionedScalar& gauge() const { return *gaugePtr_; }
    
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
            const implicitExplicit& IE
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

        // The part of the advection with Co>maxCoExp_ is solved implicitly
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
        void operator=(const MPDATA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MPDATA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
