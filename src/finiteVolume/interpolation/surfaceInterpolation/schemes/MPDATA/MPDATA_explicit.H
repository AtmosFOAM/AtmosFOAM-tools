/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    Foam::MPDATA_explicit

Description
    First apply explicit upwind interpolation to tracer T.
    Then calculate Tf = linearInterpolate(T)
    gradT = linearInterpolate(fvc::grad(T))
    Then calculate the anti-diffusive flux
    anteD = 0.5/Tf*mag(phi)*snGrad(T)/mesh.deltaCoeffs()
              - phi*dt* (Uf & gradT)
    and apply upwinding again.

SourceFiles
    MPDATA_explicit.C

\*---------------------------------------------------------------------------*/

#ifndef MPDATA_explicit_H
#define MPDATA_explicit_H

#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MPDATA_explicit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MPDATA_explicit
:
    public upwind<Type>
{
    // Private Data

        //- The MPDATA_explicit gauge value
        const scalar gauge_;

        //- Reference to the surface flux used to choose upwind direction
        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MPDATA_explicit(const MPDATA_explicit&);

        //- Disallow default bitwise assignment
        void operator=(const MPDATA_explicit&);


public:

    //- Runtime type information
    TypeName("MPDATA_explicit");


    // Constructors

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        MPDATA_explicit(const fvMesh& mesh, Istream& is)
        :
            upwind<Type>
            (
                mesh,
                mesh.lookupObject<surfaceScalarField>(word(is))
            ),
            gauge_(readScalar(is)),
            faceFlux_(upwind<Type>::faceFlux_)
        {}


        //- Construct from mesh, faceFlux and Istream
        MPDATA_explicit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux),
            gauge_(readScalar(is)),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
