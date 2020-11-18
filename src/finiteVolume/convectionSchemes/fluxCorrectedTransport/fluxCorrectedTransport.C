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

#include "fluxCorrectedTransport.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvcSurfaceIntegrateInOut.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
fluxCorrectedTransport<Type>::fluxCorrection
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& lowFlux,
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, fvPatchField, volMesh>& vfT
) const
{
    const fvMesh& mesh = this->mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    
    // Initialise the high order correction
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            interpHighScheme().interpolate(faceFlux, vf) - lowFlux
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();
    
    // High order anti-diffusive fluxes (using Zalesak FCT notation, A)
    GeometricField<Type, fvsPatchField, surfaceMesh> A
        = faceFlux*mesh.time().deltaT()*sfCorr;

    // Fluxes into and out of each cell
    volScalarField Pp = fvc::surfaceIntegrateIn(A);
    volScalarField Pm = fvc::surfaceIntegrateOut(A);
    
    // Amount each cell can rise or fall by
    volScalarField Qp = -vfT;
    volScalarField Qm = vfT;
    if (Tmax_ > Tmin_) // Assuming fixed given max and min
    {
        Qp += dimensionedScalar("", vf.dimensions(), Tmax_);
        Qm -= dimensionedScalar("", vf.dimensions(), Tmin_);
    }
    else    // Based on max and min in evolving field
    {
        volScalarField Tmin = min
        (
            vfT,
            vf.oldTime()
        );
        volScalarField Tmax = max
        (
            vfT,
            vf.oldTime()
        );
        
        for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            label own = owner[faceI];
            label nei = neighbour[faceI];

            Tmax[nei] = max(Tmax[nei], vf.oldTime()[own]);
            Tmin[own] = min(Tmin[own], vf.oldTime()[nei]);
            Tmin[nei] = min(Tmin[nei], vf.oldTime()[own]);
            Tmax[own] = max(Tmax[own], vf.oldTime()[nei]);
            if (faceFlux[faceI] >  0)
            {
                Tmax[nei] = max(Tmax[nei], max(vfT[own], vf.oldTime()[own]));
                Tmin[nei] = min(Tmin[nei], min(vfT[own], vf.oldTime()[own]));
            }
            else
            {
                Tmax[own] = max(Tmax[own], max(vfT[nei], vf.oldTime()[nei]));
                Tmin[own] = min(Tmin[own], min(vfT[nei], vf.oldTime()[nei]));
            }
        }

        Qp += Tmax;
        Qm -= Tmin;
    }
//    // Reduce Qp and Qm a bit for safety
//    Qp *= 1-SMALL;
//    Qm *= 1-SMALL;
    
    // Ratios
    volScalarField Rp = min(1., Qp/max(Pp, SMALL));
    volScalarField Rm = min(1., Qm/max(Pm, SMALL));
    forAll(Rp, cellI)
    {
        if (mag(Pp[cellI]) < VSMALL) Rp[cellI] = 0;
        if (mag(Pm[cellI]) < VSMALL) Rm[cellI] = 0;
    }
    
    // Limit the fluxes
    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];
        
        // Do not correct if the corrections are down gradient
        if (A[faceI] * (vfT[nei] - vfT[own]) < 0)
        {
            sfCorr[faceI] = 0;
        }
        else if (A[faceI] > 0)
        {
            sfCorr[faceI] *= min(Rp[nei], Rm[own]);
        }
        else
        {
            sfCorr[faceI] *= min(Rp[own], Rm[nei]);
        }
    }
    // Zero correction at boundaries
    for(label ib = 0; ib < sfCorr.boundaryField().size(); ib++)
    {
        sfCorr.boundaryFieldRef()[ib] *= 0;
    }
    
    return tsfCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
