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
    Foam::geostrophic_p_rghFvPatchScalarField

Description
    Sets pressure gradient to be 2 rho (Omega ^ Ug) at the boundary
    where Ug is the geostrophic velocity

SourceFiles
    geostrophic_p_rghFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef geostrophic_p_rghFvPatchScalarFields_H
#define geostrophic_p_rghFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class geostrophic_p_rghFvPatch Declaration
\*---------------------------------------------------------------------------*/

class geostrophic_p_rghFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data
    
        //- The geostrophic velocity
        vector Ug_;
        
        //- The rotation rate of the geometry
        vector Omega_;
    
public:

    //- Runtime type information
    TypeName("geostrophic_p_rgh");


    // Constructors

        //- Construct from patch and internal field
        geostrophic_p_rghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        geostrophic_p_rghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  geostrophic_p_rghFvPatchScalarField onto a new patch
        geostrophic_p_rghFvPatchScalarField
        (
            const geostrophic_p_rghFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        geostrophic_p_rghFvPatchScalarField
        (
            const geostrophic_p_rghFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new geostrophic_p_rghFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        geostrophic_p_rghFvPatchScalarField
        (
            const geostrophic_p_rghFvPatchScalarField&,
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
                new geostrophic_p_rghFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            // The geostrophic velocity
            const vector& Ug() const { return Ug_;}
            vector& Ug() { return Ug_;}

            // The rotation of the geometry
            const vector& Omega() const { return Omega_;}
            vector& Omega() { return Omega_;}

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
