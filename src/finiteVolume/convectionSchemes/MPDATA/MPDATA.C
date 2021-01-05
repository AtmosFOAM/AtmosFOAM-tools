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

#include "MPDATA.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"
#include "fvc.H"
#include "uncorrectedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
MPDATA<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Just use the low order interpolate for now
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            upwindConvect().interpolate(faceFlux, vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
MPDATA<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
MPDATA<Type>::fvmDiv
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
    
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    //const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Create temporary field to advect
    GeometricField<Type, fvPatchField, volMesh> vfT = vf;
    
    // Create and solve the advection equation for the low order scheme
    EulerDdtScheme<Type> backwardEuler(mesh);
    fvMatrix<Type> fvmT
    (
        backwardEuler.fvmDdt(vfT)
      + upwindConvect().fvmDiv(faceFlux, vfT)
    );
    fvmT.solve();
    
    // Put all the divergence in the RHS
    fvm += upwindConvect().fvcDiv(faceFlux, vfT);

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
MPDATA<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Calculate necessary additional fields for the correction
    // The full velocity field from the flux
    surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(faceFlux));
    Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

    // The volume field interpolated onto faces
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf = linearInterpolate(vf);

    // The surface normal gradient and the guass gradient interpolated to faces
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    GeometricField<Type, fvsPatchField, surfaceMesh> snGradT
         = snGrad.snGrad(vf);

    fv::gaussGrad<Type> grad(mesh);
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    > gradT = linearInterpolate
    (
        grad.gradf(Tf, "gradf")
    );

    // Gauge Stabilisation
    dimensionedScalar gauge("gauge", Tf.dimensions(), max(gauge_, SMALL));
    Tf += gauge;

    // Ante-diffusive flux
    surfaceScalarField anteD = 0.5/Tf*
    (
        mag(faceFlux)*snGradT/rdelta
      - faceFlux*dt*(Uf & gradT)
    );

    // Limit the anti-diffusive velocity so that Courant<0.5
    surfaceScalarField CoLim = 0.25*mesh.magSf()/rdelta/dt;
    anteD = min(anteD, CoLim);
    anteD = max(anteD, -CoLim);
    
    // Ante-diffusive velocity recontruct and then interpolate to smooth
    volVectorField V("anteDV", fvc::reconstruct(anteD));
    anteD = linearInterpolate(V) & mesh.Sf();

    // Write out the ante-diffusive velocity if needed
    if (mesh.time().writeTime())
    {
        V.write();
    }
    
    // Advect the volume field first with upwind and then the correction
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv(faceFlux, vf)
    );

    GeometricField<Type, fvPatchField, volMesh> T = vf - dt*tConvection();
    
    upwind<Type> upwindAndteD(mesh, anteD);
    tConvection.ref() += fvc::div(anteD*upwindAndteD.interpolate(T));

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
