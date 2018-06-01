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

#include "hydrostatic_p_rghFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrostatic_p_rghFvPatchScalarField::
hydrostatic_p_rghFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


hydrostatic_p_rghFvPatchScalarField::
hydrostatic_p_rghFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


hydrostatic_p_rghFvPatchScalarField::
hydrostatic_p_rghFvPatchScalarField
(
    const hydrostatic_p_rghFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


hydrostatic_p_rghFvPatchScalarField::
hydrostatic_p_rghFvPatchScalarField
(
    const hydrostatic_p_rghFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf)
{}


hydrostatic_p_rghFvPatchScalarField::
hydrostatic_p_rghFvPatchScalarField
(
    const hydrostatic_p_rghFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostatic_p_rghFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rho
         = patch().lookupPatchField<volScalarField, scalar>("rho");
    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    gradient() = -(g.value() & patch().Cf())*rho.snGrad();
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void hydrostatic_p_rghFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    hydrostatic_p_rghFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
