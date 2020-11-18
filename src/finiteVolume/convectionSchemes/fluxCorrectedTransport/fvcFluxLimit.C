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

#include "fvcFluxLimit.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fluxLimit
(
    surfaceScalarField& fluxCorr,
    const volScalarField& Tint,
    const volScalarField& Told
)
{
    const fvMesh& mesh = fluxCorr.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    
    // Fluxes into and out of each cell
    volScalarField Pp = fvc::surfaceIntegrateIn(fluxCorr);
    volScalarField Pm = fvc::surfaceIntegrateOut(fluxCorr);
    
    // Local extrema
    volScalarField Tmin = min(Tint, Told);
    volScalarField Tmax = max(Tint, Told);
        
    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        Tmax[nei] = max(Tmax[nei], max(Tint[own], Told[own]));
        Tmax[own] = max(Tmax[own], max(Tint[nei], Told[nei]));
        Tmin[nei] = min(Tmin[nei], min(Tint[own], Told[own]));
        Tmin[own] = min(Tmin[own], min(Tint[nei], Told[nei]));
    }

    // Amount each cell can rise or fall by
    volScalarField Qp = Tmax - Tint;
    volScalarField Qm = Tint - Tmin;

    // Ratios
    dimensionedScalar Psmall("Psmall", Pp.dimensions(), SMALL);
    volScalarField Rp = min(1., Qp/max(Pp, Psmall));
    volScalarField Rm = min(1., Qm/max(Pm, Psmall));
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
        if (fluxCorr[faceI] * (Tint[nei] - Tint[own]) < 0)
        {
            fluxCorr[faceI] = 0;
        }
        else if (fluxCorr[faceI] > 0)
        {
            fluxCorr[faceI] *= min(Rp[nei], Rm[own]);
        }
        else
        {
            fluxCorr[faceI] *= min(Rp[own], Rm[nei]);
        }
    }
    // Zero correction at boundaries
    for(label ib = 0; ib < fluxCorr.boundaryField().size(); ib++)
    {
        fluxCorr.boundaryFieldRef()[ib] *= 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
