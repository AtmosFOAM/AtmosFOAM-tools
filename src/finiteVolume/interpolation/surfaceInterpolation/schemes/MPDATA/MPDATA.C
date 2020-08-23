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

#include "MPDATA.H"
#include "fvc.H"
#include "uncorrectedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::MPDATA<Type>::correction
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
                "MPDATA::correction(" + vf.name() + ')',
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
         = linearInterpolate(vf.oldTime());

    upwind<Type> upInterp(mesh, faceFlux_);
    GeometricField<Type, fvsPatchField, surfaceMesh> Tfup
         = upInterp.interpolate(vf.oldTime());

    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    GeometricField<Type, fvsPatchField, surfaceMesh> snGradT
         = snGrad.snGrad(vf.oldTime());

    GeometricField<Type, fvPatchField, volMesh> Tup
         = vf.oldTime() - dt*fvc::surfaceIntegrate(faceFlux_*Tfup);

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
    //Tf = max(Tf, dimensionedScalar("", Tf.dimensions(), SMALL));
    //Tf = stabilise(Tf, dimensionedScalar("", Tf.dimensions(), SMALL));
    Tf += SMALL*max(Tf);

    // Ante-diffusive flux
    surfaceScalarField anteD = 0.5/Tf*
    (
        mag(faceFlux_)*snGradT/rdelta
      - faceFlux_*dt*(Uf & gradT)
    );

    // Limit the anti-diffusive velocity so that Courant<1
    anteD = min(anteD, 0.5*mesh.magSf()/rdelta/dt);
    anteD = max(anteD, -0.5*mesh.magSf()/rdelta/dt);

    // Correction
    upwind<Type> upwindAndteD(mesh, anteD);
    GeometricField<Type, fvsPatchField, surfaceMesh> sfCorrFlux
         = upwindAndteD.interpolate(Tup);
    sfCorr = sfCorrFlux*anteD/stabilise
    (
        faceFlux_,
        dimensionedScalar("", faceFlux_.dimensions(), SMALL)
    );

    return tsfCorr;
}


namespace Foam
{
    //makeSurfaceInterpolationScheme(MPDATA)
    makeSurfaceInterpolationTypeScheme(MPDATA, scalar)
    //makeSurfaceInterpolationTypeScheme(MPDATA, vector)
}

// ************************************************************************* //
