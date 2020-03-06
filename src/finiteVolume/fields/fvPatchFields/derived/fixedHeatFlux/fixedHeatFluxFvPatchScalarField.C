/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fixedHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedHeatFluxFvPatchScalarField::
fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatFlux_(p.size(), scalar(0)),
    densityName_("densityNameMustBeSet"),
    diffusivityName_("diffusivityNameMustBeSet")
{}


fixedHeatFluxFvPatchScalarField::
fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    heatFlux_("heatFlux", dict, p.size()),
    densityName_(dict.lookup("densityName")),
    diffusivityName_(dict.lookup("diffusivityName"))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        //gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
//        gradient() = 0.0;
    }
}


fixedHeatFluxFvPatchScalarField::
fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatFlux_(p.size(), scalar(0)),
    densityName_("densityNameMustBeSet"),
    diffusivityName_("diffusivityNameMustBeSet")
{}


fixedHeatFluxFvPatchScalarField::
fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    heatFlux_(wbppsf.heatFlux_),
    densityName_(wbppsf.densityName_),
    diffusivityName_(wbppsf.diffusivityName_)
{}


fixedHeatFluxFvPatchScalarField::
fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    heatFlux_(wbppsf.heatFlux_),
    densityName_(wbppsf.densityName_),
    diffusivityName_(wbppsf.diffusivityName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rho
         = patch().lookupPatchField<volScalarField, scalar>(densityName_);
    const uniformDimensionedScalarField& alpha =
        db().lookupObject<uniformDimensionedScalarField>(diffusivityName_);

    gradient() = heatFlux_/(rho*alpha.value());
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
    writeEntry(os, "heatFlux", heatFlux_);
    writeEntry(os, "densityName", densityName_);
    writeEntry(os, "diffusivityName", diffusivityName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedHeatFluxFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
