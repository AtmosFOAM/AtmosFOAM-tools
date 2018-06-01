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

Class
    Foam::hydrostatic_p_rghFvPatchScalarField

Description
    Sets pressure gradient of p_rgh to be (g&mesh.Cf()) grad(rho)
    for hydrostatic balance with no acceleration

SourceFiles
    hydrostatic_p_rghFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostatic_p_rghFvPatchScalarFields_H
#define hydrostatic_p_rghFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class hydrostatic_p_rghFvPatch Declaration
\*---------------------------------------------------------------------------*/

class hydrostatic_p_rghFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("hydrostatic_p_rgh");


    // Constructors

        //- Construct from patch and internal field
        hydrostatic_p_rghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrostatic_p_rghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hydrostatic_p_rghFvPatchScalarField onto a new patch
        hydrostatic_p_rghFvPatchScalarField
        (
            const hydrostatic_p_rghFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        hydrostatic_p_rghFvPatchScalarField
        (
            const hydrostatic_p_rghFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new hydrostatic_p_rghFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        hydrostatic_p_rghFvPatchScalarField
        (
            const hydrostatic_p_rghFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new hydrostatic_p_rghFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
