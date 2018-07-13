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
    Foam::Partitioned

Description
    A list of PtrList of things, one for each partion of the atmosphere. 
    It is also an MeshObject so that things can be looked up from the database.
    It holds a wordList for the partition names.

SourceFiles
    Partitioned.C

\*---------------------------------------------------------------------------*/

#ifndef Partitioned_H
#define Partitioned_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Partitioned Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Partitioned
:
    public MeshObject<fvMesh, MoveableMeshObject, Partitioned<Type> >,
    PtrList<Type>
{
    // Private data

        //- The names of the parts
        const wordList partNames_;

    // Dissallow copy construct
        Partitioned(const Partitioned&);


public:

    // Constructors

        //- Construct empty
        Partitioned(const fvMesh& mesh);

        //- Construct from components
        Partitioned
        (
            const fvMesh& mesh,
            const wordList& partNames__
        );

    //- Destructor
    ~Partitioned();

    // Access functions
    
        //- The names of the parts
        const wordList& partNames() const { return partNames_; }


    // Member functions
        //- No changes needed for moving the mesh
        bool movePoints() {return true;}

    // Functions from PtrList

        label size() const
        {
            return PtrList<Type>::size();
        }
        
        const Type& operator[]
        (
            const label ip
        ) const
        {
            return PtrList<Type>::operator[](ip);
        }

        Type& operator[]
        (
            const label ip
        )
        {
            return PtrList<Type>::operator[](ip);
        }
        
        autoPtr<Type> set(label ip, autoPtr<Type> T)
        {
            return PtrList<Type>::set(ip, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Partitioned.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
