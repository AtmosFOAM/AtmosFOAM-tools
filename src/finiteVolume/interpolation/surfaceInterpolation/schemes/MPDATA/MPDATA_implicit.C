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

\*---------------------------------------------------------------------------*/

#include "MPDATA_implicit.H"
#include "fvc.H"
#include "uncorrectedSnGrad.H"
#include "CourantNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::MPDATA_implicit<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "MPDATA_implicit::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Calculate necessary additional fields
    surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(faceFlux_));
    Uf += (faceFlux_ - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf
         = linearInterpolate(vf);

//    upwind<Type> upInterp(mesh, faceFlux_);
//    GeometricField<Type, fvsPatchField, surfaceMesh> Tfup
//         = upInterp.interpolate(vf);

    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    GeometricField<Type, fvsPatchField, surfaceMesh> snGradT
         = snGrad.snGrad(vf);

//    GeometricField<Type, fvPatchField, volMesh> Tup
//         = vf.oldTime() - dt*fvc::surfaceIntegrate(faceFlux_*Tfup);

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

    // Stabilisation
    Tf = max(Tf, dimensionedScalar("", Tf.dimensions(), SMALL));

    // Ante-diffusive flux
    surfaceScalarField anteD = 0.5/Tf*
    (
        mag(faceFlux_)*snGradT/rdelta
      + faceFlux_*dt*(Uf & gradT)
    );
    
    // Limit the anti-diffusive velocity so that Courant<0.5
    surfaceScalarField CoLim = 0.25*mesh.magSf()/rdelta/dt;
    anteD = min(anteD, CoLim);
    anteD = max(anteD, -CoLim);
    
    // Ante-diffusive velocity recontruct and then interpolate to smooth
    volVectorField V("anteDV", fvc::reconstruct(anteD));
    anteD = linearInterpolate(V) & mesh.Sf();

    // Calculate and write out the ante-diffusive flux and Co if needed
    if (mesh.time().writeTime())
    {
        V.write();
        volScalarField Co("anteDCo", CourantNo(anteD, dt));
        Co.write();
    }

    // Correction
    upwind<Type> upwindAndteD(mesh, anteD);
    GeometricField<Type, fvsPatchField, surfaceMesh> sfCorrFlux
         = upwindAndteD.interpolate(vf);
    sfCorr = sfCorrFlux*anteD/stabilise
    (
        faceFlux_,
        dimensionedScalar("", faceFlux_.dimensions(), SMALL)
    );

    return tsfCorr;
}


namespace Foam
{
    //makeSurfaceInterpolationScheme(MPDATA_implicit)
    makeSurfaceInterpolationTypeScheme(MPDATA_implicit, scalar)
    //makeSurfaceInterpolationTypeScheme(MPDATA_implicit, vector)
}

// ************************************************************************* //
