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

#include "splitConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
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
splitConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Initialise the surface field as the low order interpolate
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            upwindConvect().interpolate(faceFlux, vf)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& interp = tinterp.ref();
    
    // Low order advection of the vol field
    GeometricField<Type, fvPatchField, volMesh> Ttd = vf.oldTime()
        - mesh.time().deltaT()*fvc::div(faceFlux*interp);

    // High order correction for the interpolation
    interp += fluxCorrection(interp, faceFlux, vf, Ttd);

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
splitConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
splitConvectionScheme<Type>::fvmDiv
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

    // Create temporary field to advect
    GeometricField<Type, fvPatchField, volMesh> vfT(vf);
    vfT.oldTime() = vf.oldTime();
    
    // Create and solve the advection equation for the low order scheme
    EulerDdtScheme<Type> backwardEuler(this->mesh());
    fvMatrix<Type> fvmT
    (
        backwardEuler.fvmDdt(vfT)
      + upwindConvect().fvmDiv(faceFlux, vfT)
    );
    fvmT.solve();
    
    // Put all the divergence in the RHS
    GeometricField<Type, fvsPatchField, surfaceMesh> lowFlux
        = upwindConvect().interpolate(faceFlux, vfT);
    fvm.source() = fvc::surfaceIntegrate
    (
        faceFlux*(lowFlux + fluxCorrection(lowFlux, faceFlux, vf, vfT))
    );

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
splitConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
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
