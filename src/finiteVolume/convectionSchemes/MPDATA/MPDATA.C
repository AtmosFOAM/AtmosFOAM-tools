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
#include "CourantNoFunc.H"
#include "localMax.H"
#include "leastSquaresGrad.H"
#include "fvcFluxLimit.H"
#include "fvcLocalMinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void MPDATA<Type>::calculateAnteD
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

    // The volume field interpolated onto faces for the denominator
    /*volScalarField T = max
    (
        vf + dimensionedScalar("", vf.dimensions(), gauge_ + SMALL),
        8*dt*mag(fvc::div(faceFlux, vf, "MPDATA_div"))
    );*/
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf
         = fvc::interpolate(vf, "MPDATA_denom");

    // Stabilisation
    Tf += dimensionedScalar("", Tf.dimensions(), gauge_ + SMALL);

    // Gradient of T in the cell centre to cell centre direction
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    surfaceScalarField snGradT = snGrad.snGrad(vf);

    // The correction in space
    anteD() = 0.5/Tf*mag(faceFlux)*snGradT/rdelta;
    
    // Apply a correction in time if needed
    if (timeCorrector_ == "advective")
    {
        // The full velocity field from the flux and correct
        surfaceVectorField Uf = fvc::interpolate
        (
            fvc::reconstruct(faceFlux), "MPDATA_velocity"
        );
        Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        // Full face gradient of T
        fv::leastSquaresGrad<Type> grad(mesh);
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        > gradT = fvc::interpolate(grad.grad(vf), "MPDATA_gradient");
        gradT += (snGradT - (gradT & mesh.delta())*rdelta)
             * mesh.delta()*rdelta;

        // add advetive form time correction
        anteD() -= max(1-2*offCentre, scalar(0))*0.5/Tf*faceFlux*dt*(Uf & gradT);
        
        // reduce for large gradients in offCentre
        //localMax<vector> maxInterp(this->mesh());
        //surfaceVectorField gradOff(maxInterp.interpolate(fvc::grad(offCentre)));
        //anteD() /= 1 + 10*dt*mag(gradOff)*mag(Uf);
        //anteD() /= 1 + dt*mag(linearInterpolate(fvc::div(faceFlux*offCentre)));
    }
    else if(timeCorrector_ == "flux")
    {
        anteD() -= max(1-2*offCentre, scalar(0))*0.5/Tf*faceFlux*dt
         *fvc::interpolate(fvc::div(faceFlux, vf, "MPDATA_div"), "MPDATA_idiv");
    }
    
    // Smooth where offCentre>0
    surfaceVectorField V("anteDV", linearInterpolate(fvc::reconstruct(anteD())));
    surfaceScalarField imp = min(2*offCentre, scalar(1));
    imp = maxInterp.interpolate(fvc::localMax(imp));
    imp = linearInterpolate(fvc::localMax(imp));
    anteD() = imp*(V & mesh.Sf()) + (1-imp)*anteD();
    
    // Limit to obey Courant number restriction
    volScalarField CoV("CoV", 4*CourantNo(anteD(), dt));
    //surfaceScalarField Cof("Cof", 8*dt*mag(anteD())/faceVol_);
    //Cof = maxInterp.interpolate(fvc::localMax(Cof));
    //Cof = linearInterpolate(fvc::localMax(Cof));
    surfaceScalarField Cof("Cof", maxInterp.interpolate(CoV));
    anteD() /= max(scalar(1), Cof);
    /*if (mesh.time().writeTime())
    {
        Cof.write();
        CoV.write();
    }*/

    //anteD() = sign(anteD())*min(mag(anteD()), 0.125*faceVol_/dt);

/*    // Write out the ante-diffusive velocity if needed
    if (mesh.time().writeTime())
    {
        // Ante-diffusive velocity and divergence
        surfaceVectorField V
        (
            "anteDV",
            linearInterpolate(fvc::reconstruct(anteD()))
        );
        V += (anteD() - (V & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        V.write();
        
        // divergence of V
        volScalarField divAnteD("divAnteD", fvc::div(anteD()));
        divAnteD.write();
    }*/
    CoV = CourantNo(anteD(), dt);
    Info << "Ante-diffusive Courant number max: " << max(CoV).value() << endl;

}

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
    tfvm.ref() += fvcDiv(faceFlux, vf);
    
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

    //localMax<scalar> maxInterp(this->mesh());
    //volScalarField C = CourantNo(faceFlux, dt);
    surfaceScalarField Cf = 2*dt*mag(faceFlux)/faceVol_;/*max
    (
        2*dt*mag(faceFlux)/faceVol_, 
        maxInterp.interpolate(C)
    );*/
    // Smooth Cf to get smooth offCentre
    //Cf = maxInterp.interpolate(fvc::localMax(Cf));
    //Cf = linearInterpolate(fvc::localMax(Cf));
    const surfaceScalarField offCentre = offCentre_ < 0 ?
        surfaceScalarField("offCentre", max(1-1/(Cf + 0.2), scalar(0))) :
        surfaceScalarField
        (
            IOobject("offCentre", mesh.time().timeName(), mesh),
            mesh,
            dimensionedScalar("", dimless, offCentre_),
            "fixedValue"
        );
    Info << "offCentre goes from " << min(offCentre).value() << " to " 
         << max(offCentre).value() << endl;
    
    if (mesh.time().writeTime())
    {
        offCentre.write();
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
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf - dt*tConvection(),
        vf.boundaryField().types()
    );
    
    // Calculte the implicit part of the advection
    if (mag(offCentre_) > SMALL)
    {
        EulerDdtScheme<Type> backwardEuler(this->mesh());
        fvMatrix<Type> fvmT
        (
            backwardEuler.fvmDdt(T)
          + upwindConvect().fvmDiv(offCentre*faceFlux, T)
        );
        fvmT.solve();
    
        // Add the low order implicit divergence
        tConvection.ref() += upwindConvect().fvcDiv(offCentre*faceFlux, T);
    }
    
    // Calculate, apply (and update) the correction
    //if (nCorr_ > 0) calculateAnteD(faceFlux, vf, offCentre);
    if (nCorr_ > 0) calculateAnteD(faceFlux, T, offCentre);
    for(label iCorr =1; iCorr < nCorr_; iCorr++)
    {
        calculateAnteD
        (
            faceFlux,
            T - dt*anteDConvect().fvcDiv(anteD(), T),
            offCentre
        );
    }
    if (nCorr_ > 0)
    {
        dimensionedScalar gauge("gauge", T.dimensions(), gauge_);
        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_+gauge_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_+gauge_);
        T += gauge;
        anteD() *= anteDConvect().interpolate(anteD(), T);
        
        // Limit the fluxes with Zalesak FCT limiter if needed
        if (FCTlimit_)
        {
            if (FCTmin_ < FCTmax_)
            {
                fvc::fluxLimit(anteD(), T, FCTmin, FCTmax, dt);
            }
            else if (mag(offCentre_) < SMALL)
            {
                fvc::fluxLimit(anteD(), T, vf+gauge, dt);
            }
            else
            {
                fvc::fluxLimit(anteD(), T, dt);
            }
        }
    
        tConvection.ref() += fvc::div(anteD());
    }
    
    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
