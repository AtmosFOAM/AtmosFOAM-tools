/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Class
    Foam::gaussian

Description
    A class which describes a Gaussian fucntion based on its maxima, radius
    and centre. Returns the Gaussian as a tmp<volScalarField>

SourceFiles
    gaussian.C
    gaussianIO.C
    
\*---------------------------------------------------------------------------*/

#ifndef gaussian_H
#define gaussian_H

#include "volFieldsFwd.H"
#include "point.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
class fvMesh;

// Forward declaration of friend functions and operators
class gaussian;
Istream& operator>>(Istream&, gaussian&);
Ostream& operator<<(Ostream&, const gaussian&);

/*---------------------------------------------------------------------------*\
                         Class gaussian Declaration
\*---------------------------------------------------------------------------*/

class gaussian
{
    // Private data

        //- Centre of the Gaussian
        dimensioned<point> centre_;
        
        //- Radius (standard deviation) of the Gaussian (in metres)
        dimensionedScalar radius_;
        
        //- Maximum value of the Gaussian (at the centre)
        dimensionedScalar max_;

        //- Location of the periodic minimum of the domain
        point periodMin_;

        //- Location of the periodic maximum of the domain
        point periodMax_;

public:

    // Constructors

        //- Construct null
        gaussian();

        //- Construct from components
        gaussian
        (
            const dimensioned<point>& centre__,
            const dimensionedScalar radius__,
            const dimensionedScalar max__,
            const point& periodMin__ = point::zero,
            const point& periodMax__ = point::zero
        );

        //- Construct from Istream
        gaussian
        (
            Istream&,
            const point& periodMin__ = point::zero,
            const point& periodMax__ = point::zero
        );

        //- Construct as copy
        gaussian(const gaussian&);


    //- Destructor
    ~gaussian();


    // Member Functions

        // Access (const and non-const access to all members)
        const dimensioned<point>& centre() const { return centre_;}
        dimensioned<point>& centre() { return centre_;}
        const dimensionedScalar& radius() const { return radius_;}
        dimensionedScalar& radius() { return radius_;}
        const dimensionedScalar& max() const { return max_;}
        dimensionedScalar& max() { return max_;}
        const point& periodMin() const { return periodMin_;}
        point& periodMin() { return periodMin_;}
        const point& periodMax() const { return periodMax_;}
        point& periodMax() { return periodMax_;}

        // Return a volScalarField for the Gaussian distribution
        tmp<volScalarField> field(const fvMesh& mesh) const;

    // Member Operators

        void operator=(const gaussian&);

    // IOstream Operators

        friend Istream& operator>>(Istream&, gaussian&);
        friend Ostream& operator<<(Ostream&, const gaussian&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
