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
    const implicitExplicit& IE
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Calculate necessary additional fields for the correction
    // The full velocity field from the flux and correct
    surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(faceFlux));
    Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

    // The volume field interpolated onto faces
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf = linearInterpolate(vf);

    // The volume field averaged back onto cell centres
    GeometricField<Type, fvPatchField, volMesh> Tc = fvc::average(Tf);

    // The surface normal gradient and the guass gradient interpolated to faces
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    GeometricField<Type, fvsPatchField, surfaceMesh> snGradT
         = snGrad.snGrad(vf);

    // Gauge Stabilisation
    gauge().dimensions() = Tf.dimensions();
    Tf += gauge() + dimensionedScalar("", Tf.dimensions(), VSMALL);
    Tc += gauge() + dimensionedScalar("", Tf.dimensions(), VSMALL);

    fv::gaussGrad<Type> grad(mesh);
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    > gradTbyT = linearInterpolate
    (
        grad.gradf(Tf, "gradf")/Tc
    );
    gradTbyT += (snGradT/Tf - (gradTbyT & mesh.Sf())/mesh.magSf())
         * mesh.Sf()/mesh.magSf();

    int IEalpha = IE == EXPLICIT ? -1 : 1;

    // Smooth by increasing Tf for implicit advection
    if (IE == IMPLICIT)
    {
        surfaceScalarField minTf = 2*dt/mesh.magSf()*mag
        (
            mag(faceFlux)*snGradT
          + IEalpha*faceFlux*dt*(Uf & gradTbyT)*Tf*rdelta
        );
        minTf = linearInterpolate(fvc::average(minTf));
        
        //Tf = max(Tf, minTf);
        //Tf += minTf;
        Tf = sqrt(sqr(Tf) + sqr(minTf));
    }

    // ante-diffusive flux for implicit or explicit advection
    anteD() = 0.5*
    (
        mag(faceFlux)*snGradT/(Tf*rdelta)
     + IEalpha*faceFlux*dt*(Uf & gradTbyT)
//         + IEalpha*faceFlux*dt*linearInterpolate(fvc::div(faceFlux*Tf)/Tc)
    );
    
    // Limit the anti-diffusive Courant number for implicit advection
/*    if (IE == IMPLICIT)
    {
        volScalarField Co("anteDeCo", CourantNo(anteD(), dt));
        surfaceScalarField Cof = linearInterpolate(Co);
        anteD() *= min(Cof, scalar(0.5))/(Cof + VSMALL);
    }
*/
    // Ante-diffusive velocity recontruct (and then interpolate to smooth)
    word Vname = IE == EXPLICIT? "anteDVe" : "anteDVi";
    volVectorField V(Vname, fvc::reconstruct(anteD()));
    //anteD() = linearInterpolate(V) & mesh.Sf();
    //upwind<vector> dInterp(mesh,faceFlux);
    //anteD() = dInterp.interpolate(V) & mesh.Sf();

    // Write out the ante-diffusive velocity if needed
    if (mesh.time().writeTime())
    {
        V.write();
        volScalarField Co("anteDeCo", CourantNo(anteD(), dt));
        Co.write();
    }
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
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Separate faceFlux into explicit and implicitly solved parts
    surfaceScalarField fluxSmall = 0*faceFlux;
    surfaceScalarField fluxBig = faceFlux;
    if (maxCoExp_>0)
    {
        fluxSmall = faceFlux;
        fluxBig *= 0;
    
        surfaceScalarField fluxLimit = 0.5*maxCoExp_*mesh.magSf()/rdelta/dt;
        forAll(faceFlux, faceI)
        {
            if (mag(faceFlux[faceI]) > fluxLimit[faceI])
            {
                fluxSmall[faceI] = 0;
                fluxBig[faceI] = faceFlux[faceI];
            }
        }
    }

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
    
    // Create temporary field to advect and the temporary divergence field
    GeometricField<Type, fvPatchField, volMesh> T = vf;
    
    // The temporary divergence field is initialised as the explicit divergence
    GeometricField<Type, fvPatchField, volMesh> div = fvcDiv(fluxSmall, vf);
    
    // First apply the explicit part
    if (maxCoExp_>0)
    {
        // Add the explicit MPDATA correction to the RHS
        fvm += div;
        T -= dt*div;
    }
    
    // Next calculate and apply the correction for the implicit part
    if (implicitCorrection_ > 0)
    {
        // The advection equation for the imlicit high order correction
        calculateAnteD(fluxBig, vf, IMPLICIT);

        // Add the implicit MPDATA correction to the RHS
        div = anteDConvect().fvcDiv(implicitCorrection_*anteD(), T+gauge());
        fvm += div;
        T -= dt*div;
    }

    // implicit low order part to update T from T.oldTime()
    T.oldTime() = T;
    EulerDdtScheme<Type> backwardEuler(this->mesh());
    fvMatrix<Type> fvmT
    (
        backwardEuler.fvmDdt(T)
      + upwindConvect().fvmDiv(fluxBig, T)
    );
    fvmT.solve();
    // Put the low order implicit divergence on the RHS of the matrix
    fvm += upwindConvect().fvcDiv(fluxBig, T);
    
    // Re-calculate and apply the correction and limit the Courant number
    if (implicitCorrection_ > 0)
    {
        // Save the old flux
        surfaceScalarField oldCorr = anteD();
        // Re-calculate the correction
        calculateAnteD(fluxBig, T, IMPLICIT);
        anteD() -= oldCorr;
        
        surfaceScalarField c = linearInterpolate(CourantNo(fluxBig, dt));
        //anteD() *= min(c, scalar(0.5))/(c + VSMALL);
        c = sqr(min(1/(c+SMALL), dimensionedScalar("",dimless, scalar(1))));
        //c = 0.5*(1 + tanh((0.75-c)*8));
        //anteD() = c*anteD();

        fvm += anteDConvect().fvcDiv(implicitCorrection_*anteD(), T+gauge());
    }
    
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
    const dimensionedScalar& dt = this->mesh().time().deltaT();

    // Initialise the divergence to be the first-order upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv(faceFlux, vf)
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    // Store the original and upwind advected tracer
    GeometricField<Type, fvPatchField, volMesh> T = vf - dt*tConvection();
    
//    // Apply the correction (once)
//    calculateAnteD(faceFlux, vf, EXPLICIT);
//    tConvection.ref() += anteDConvect().fvcDiv(anteD(), T+gauge());
    
    // Apply multiple corrections
    T.oldTime() = vf;
    for (label iCorr = 0; iCorr < nCorr_; iCorr++)
    {
        calculateAnteD(faceFlux, T.oldTime(), EXPLICIT);
        tConvection.ref() = anteDConvect().fvcDiv(anteD(), T+gauge());
        T.oldTime() = T - dt*tConvection.ref();
    
        tConvection.ref() == (vf - T.oldTime())/dt;
    }
    
    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
