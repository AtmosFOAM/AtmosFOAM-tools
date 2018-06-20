/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::PartitionedField

Description
    A list of GeometricFields, one for each partion of the atmosphere. 
    Holds a reference to an optional partionedField sigma which is the 
    fraction of the total so that 
    sum = sum(sigma*field)
    and 
    mean = sum(sigma*field)/sum(sigma)
    or if sigma does not exist then
    sum = sum(field)

SourceFiles
    PartitionedField.C

\*---------------------------------------------------------------------------*/

#ifndef PartitionedField_H
#define PartitionedField_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PartitionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class PartitionedField
:
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Private data

        //- The base name for all of the partitions
        const word baseName_;

        //- The names of the parts
        const wordList partNames_;

        //- The sum over all partitions (weighted by sigma)
        GeometricField<Type, PatchField, GeoMesh> sum_;

        //- The mean over all partitions (weighted by sigma)
        GeometricField<Type, PatchField, GeoMesh> mean_;

        //- Switch - does the PartitionedField hold a reference to sigma
        const Switch needsSigma_;

        //- The fractions in each partition
        const PartitionedField<scalar, PatchField, GeoMesh>* sigmaPtr_;

        //- Optional storage for the rates of change
        mutable PtrList<GeometricField<Type, PatchField, GeoMesh> >* ddtPtr_;

    // Dissallow copy construct
        //PartitionedField(const PartitionedField&);


public:

    // Public typedefs

        //- Type of mesh on which this PartitionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;
        
    // Constructors

        //- Read in based on names
        PartitionedField
        (
            const word& baseName__,
            const wordList& partNames__,
            const Mesh& mesh,
            const word& timeName
        );

        //- Read in based on names and sigma
        PartitionedField
        (
            const word& baseName__,
            const wordList& partNames__,
            const Mesh& mesh,
            const word& timeName, 
            const PartitionedField<scalar, PatchField, GeoMesh>& sigma__
        );

        //- Initialise all parts to "field"
        PartitionedField
        (
            const word& baseName__,
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field
        );

        //- Initialise all parts to "field", with sigma
        PartitionedField
        (
            const word& baseName__,
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field,
            const PartitionedField<scalar, PatchField, GeoMesh>& sigma__
        );


    //- Destructor
    ~PartitionedField();

    // Access functions
    
        //- The base name for all of the partitions
        const word& baseName() const { return baseName_; }

        //- The names of the parts
        const wordList& partNames() const { return partNames_; }

        //- The sum over all partitions (weighted by sigma)
        const GeometricField<Type, PatchField, GeoMesh>& sum() const
            { return sum_; }

        //- The mean over all partitions (weighted by sigma)
        const GeometricField<Type, PatchField, GeoMesh>& mean() const
            { return mean_; }

        //- Switch - does the PartitionedField hold a reference to sigma
        const Switch& needsSigma() const { return needsSigma_;}

        //- Return sigma
        const PartitionedField<scalar, PatchField, GeoMesh>& sigma() const
            { return *sigmaPtr_; }
            
        //- Rate of change - const access
        const PtrList<GeometricField<Type, PatchField, GeoMesh>>& ddt() const
        {
            if (!ddtPtr_)
            {
                FatalErrorIn("PartitionedField::ddt()")
                    << "ddtPtr_ not set by storeTime()"
                    << abort(FatalError);
            }
            return *ddtPtr_;
        }

        //- Rate of change - non const access
        PtrList<GeometricField<Type, PatchField, GeoMesh>>& ddt()
        {
            if (!ddtPtr_)
            {
                FatalErrorIn("PartitionedField::ddt()")
                    << "ddtPtr_ not set by storeTime()"
                    << abort(FatalError);
            }
            return *ddtPtr_;
        }

    // Member Functions

        //- Update sum and mean and return sum
        const GeometricField<Type, PatchField, GeoMesh>& updateSum();
        
        //- Update sum and mean and return mean
        const GeometricField<Type, PatchField, GeoMesh>& updateMean();
        
        //- Multiply the field by sigma
        PartitionedField<Type, PatchField, GeoMesh> timesSigma() const;

        //- Divide the field by given sigma
        PartitionedField<Type, PatchField, GeoMesh> divideBy
        (
            const PartitionedField<scalar, PatchField, GeoMesh>& sigma__,
            const word newBaseName = ""
        ) const;
        
        //- Store the old time values and the rates of change
        void storeTime();

        // Output
        void write();

        // Input (re-read all fields)
        void readUpdate();

    // Member operators
    
        void operator=(const PartitionedField<Type, PatchField, GeoMesh>&);
        void operator+=(const PartitionedField<Type, PatchField, GeoMesh>&);
        void operator-=(const PartitionedField<Type, PatchField, GeoMesh>&);

    // Functions from PtrList

        label size() const
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::size();
        }
        
        const GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label ip
        ) const
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](ip);
        }

        GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label ip
        )
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](ip);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PartitionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
