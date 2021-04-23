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

#include "MPDATA_CN.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"
#include "fvc.H"
#include "uncorrectedSnGrad.H"
#include "CourantNoFunc.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void MPDATA_CN<Type>::calculateAnteD
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& offCentre
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();
    localMax<scalar> maxInterp(this->mesh());

    // Calculate necessary additional fields for the correction

    // The full velocity field from the flux and correct
    surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(faceFlux));
    Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

    // The surface normal gradient and the guass gradient interpolated to faces
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    GeometricField<Type, fvsPatchField, surfaceMesh> snGradT
         = snGrad.snGrad(vf);

    // The volume field interpolated onto faces
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf = linearInterpolate(vf);

    // Stabilisation
    Tf += dimensionedScalar("", Tf.dimensions(), SMALL);

    // The volume field averaged back onto cell centres
    //GeometricField<Type, fvPatchField, volMesh> Tc = fvc::average(Tf);

    // Full face gradient of T /T
    fv::gaussGrad<Type> grad(mesh);
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    > gradT = linearInterpolate
    (
        grad.gradf(Tf, "gradf")
    );
    gradT += (snGradT - (gradT & mesh.Sf())/mesh.magSf())
         * mesh.Sf()/mesh.magSf();

    // Smooth by increasing Tf for implicit advection
    surfaceScalarField minTf = 4*(1+offCentre)*dt/mesh.magSf()*
    (
        mag(faceFlux)*snGradT
      - max(1-2*offCentre, scalar(0))*faceFlux*dt*(Uf & gradT)*rdelta
    );
    //minTf = maxInterp.interpolate(2*fvc::average(minTf));
    Tf = max(Tf, mag(minTf));
    //Tf = sqrt(sqr(Tf) + sqr(minTf));
    //Tf += minTf;

    // ante-diffusive flux for implicit or explicit advection
    /*anteD() = 0.5/Tf*
    (
        mag(faceFlux)*snGradT/rdelta
      - max(1-2*offCentre, scalar(0))*faceFlux*dt*(Uf & gradT)
    );*/
    anteD() = 0.5/(4*dt)*mesh.magSf()/rdelta*minTf/Tf;

    // Ante-diffusive velocity recontruct (and then interpolate to smooth)
    volVectorField V("anteDV", fvc::reconstruct(anteD()));

    // Smooth for large Courant numbers
    if (max(offCentre).value() > SMALL)
    {
        surfaceScalarField needsSmoothing = min(2*offCentre, scalar(1));
        anteD() = (1-needsSmoothing)*anteD()
                + needsSmoothing*(1-needsSmoothing)*(linearInterpolate(V) & mesh.Sf());
    }

    volScalarField Co("anteDeCo", CourantNo(anteD(), dt));
    Info << "Anti diffusive Courant number max = " << max(Co).value() << endl;

    /*// Limit the anti-diffusive flux so that Co < 0.5
    const surfaceScalarField Cof = maxInterp.interpolate(Co);
    anteD() *= min(Cof, scalar(0.45))/max(Cof, SMALL);
*/
    // Write out the ante-diffusive velocity if needed
    if (mesh.time().writeTime())
    {
        V.write();
    }
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
MPDATA_CN<Type>::interpolate
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
MPDATA_CN<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
MPDATA_CN<Type>::fvmDiv
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
    tfvm.ref() += fvcDiv(faceFlux, vf);
    
    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
MPDATA_CN<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const dimensionedScalar& dt = this->mesh().time().deltaT();
    const volScalarField offCentreC
    (
        "offCentreC",
        max(1-1/CourantNo(faceFlux, dt), scalar(0))
    );
    localMax<scalar> maxInterp(this->mesh());
    const surfaceScalarField offCentre = maxInterp.interpolate(offCentreC);

    // Write out the off centering if needed
    if (this->mesh().time().writeTime())
    {
        offCentreC.write();
    }

    // Initialise the divergence to be the first-order upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv((1-offCentre)*faceFlux, vf)
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );
    
    // Create temporary field to advect and the temporary divergence field
    GeometricField<Type, fvPatchField, volMesh> T = vf;

    // Calculte the implicit part of the advection
    EulerDdtScheme<Type> backwardEuler(this->mesh());
    fvMatrix<Type> fvmT
    (
        backwardEuler.fvmDdt(T)
      + tConvection()
      + upwindConvect().fvmDiv(offCentre*faceFlux, T)
    );
    fvmT.solve();
    
    // Add the low order implicit divergence
    tConvection.ref() += upwindConvect().fvcDiv(offCentre*faceFlux, T);

    // Calculate, apply (and update) the correction
    calculateAnteD(faceFlux, vf, offCentre);
    for(label iCorr =1; iCorr < nCorr_; iCorr++)
    {
        calculateAnteD
        (
            faceFlux,
            T - dt*anteDConvect().fvcDiv(anteD(), T),
            offCentre
        );
    }
    tConvection.ref() += anteDConvect().fvcDiv(anteD(), T);
    
    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
