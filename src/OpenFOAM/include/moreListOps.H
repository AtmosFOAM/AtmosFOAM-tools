/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

InNamspace
    Foam

Description
    More functions to operate on Lists.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef moreListOps_H
#define moreListOps_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Pre add to all members of the list
template<class T>
List<T> add(const T& a, const List<T>& list)
{
    List<T> rVal(list.size());
    
    forAll(list, i)
    {
        rVal[i] = a + list[i];
    }
    
    return rVal;
}


//- Post add to all members of the list
template<class T>
List<T> add(const List<T>& list, const T& a)
{
    List<T> rVal(list.size());
    
    forAll(list, i)
    {
        rVal[i] = list[i] + a;
    }
    
    return rVal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

