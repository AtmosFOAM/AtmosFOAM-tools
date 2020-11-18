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

\*---------------------------------------------------------------------------*/

#include "CrankNicolsonConvection.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
CrankNicolsonConvection<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Initialise the surface field as the low order interpolate
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            highConvection().interpolate(faceFlux, vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
CrankNicolsonConvection<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
CrankNicolsonConvection<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Matrix to be returned will only have a source term
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();
    
    // Reference to the time step
    const dimensionedScalar& dt = vf.mesh().time().deltaT();

    surfaceScalarField betaFlux = (1-offCenter_)*faceFlux.oldTime();
    surfaceScalarField alphaFlux = offCenter_*faceFlux;

    // Create temporary field to advect
    GeometricField<Type, fvPatchField, volMesh> vfT = vf;
    GeometricField<Type, fvPatchField, volMesh> oldDiv
         = highConvection().fvcDiv(betaFlux, vf.oldTime());
    vfT.oldTime() = vf.oldTime() - dt*oldDiv;
    
    // Create and solve the advection equation for the low order scheme
    EulerDdtScheme<Type> backwardEuler(this->mesh());
    fvMatrix<Type> fvmT
    (
        backwardEuler.fvmDdt(vfT)
      + upwindConvect().fvmDiv(alphaFlux, vfT)
    );
    for(label icorr = 0; icorr < nCorr_; icorr++)
    {
        solve
        (
            fvmT == -highConvection().fvcDiv(alphaFlux, vfT)
                    + upwindConvect().fvcDiv(alphaFlux, vfT)
        );
    }
    
    // Put all the divergence in the RHS
    fvm += oldDiv + highConvection().fvcDiv(alphaFlux, vfT);

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
CrankNicolsonConvection<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the time step
    const dimensionedScalar& dt = vf.mesh().time().deltaT();

    surfaceScalarField betaFlux = (1-offCenter_)*faceFlux.oldTime();
    surfaceScalarField alphaFlux = offCenter_*faceFlux;

    GeometricField<Type, fvPatchField, volMesh> vfT = vf;
    GeometricField<Type, fvPatchField, volMesh> oldDiv
         = highConvection().fvcDiv(betaFlux, vf.oldTime());
    vfT.oldTime() = vf.oldTime() - dt*oldDiv;
    
    for(label icorr = 0; icorr < nCorr_; icorr++)
    {
        vfT = vfT.oldTime() - dt*highConvection().fvcDiv(alphaFlux, vfT);
    }
    
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        oldDiv
      + highConvection().fvcDiv(alphaFlux, vfT)
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
