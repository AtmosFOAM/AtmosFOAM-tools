/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "primitiveMesh.H"
#include "cell.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primitiveMesh::pointInCellBB
(
    const point& p,
    label celli,
    scalar inflationFraction
) const
{
    boundBox bb
    (
        cells()[celli].points
        (
            faces(),
            points()
        ),
        false
    );

    if (inflationFraction > SMALL)
    {
        vector inflation = inflationFraction*vector::one*mag(bb.span());
        bb = boundBox(bb.min() - inflation, bb.max() + inflation);
    }

    return bb.contains(p);
}


bool Foam::primitiveMesh::pointInCell(const point& p, label celli) const
{
    const labelList& f = cells()[celli];
    const labelList& owner = this->faceOwner();
    const vectorField& cf = faceCentres();
    const vectorField& Sf = faceAreas();

    scalar magSqrp = magSqr(p);
    scalar smallProj = 0.1*magSqrp;
    scalar SMALLProj = 1e-6*magSqrp;
    //scalar VSMALLProj = SMALL*magSqrp;

    bool inCell = true;

    for(label facei = 0; facei < f.size() && inCell; facei++)
    {
        label nFace = f[facei];
        point Cf = cf[nFace];
        vector proj = p - cf[nFace];
        vector normal = Sf[nFace];
        if (owner[nFace] != celli)
        {
            normal = -normal;
        }
        scalar Cfn = Cf & normal;
        if (mag(Cfn) < smallProj)
        {
            inCell = inCell && ((normal & proj) <= SMALLProj);
        }
        else
        {
            scalar magSqrDiff = magSqr(Cf) - magSqr(p);
            if (mag(magSqrDiff) > SMALLProj)
            {
                inCell = inCell && (sign(magSqrDiff) == sign(Cfn));
            }
        }
    }

    return inCell;
}


Foam::label Foam::primitiveMesh::findNearestCell(const point& location) const
{
    const vectorField& centres = cellCentres();

    label nearestCelli = 0;
    scalar minProximity = magSqr(centres[0] - location);

    for (label celli = 1; celli < centres.size(); celli++)
    {
        scalar proximity = magSqr(centres[celli] - location);

        if (proximity < minProximity)
        {
            nearestCelli = celli;
            minProximity = proximity;
        }
    }

    return nearestCelli;
}


Foam::label Foam::primitiveMesh::findCell(const point& location) const
{
    if (nCells() == 0)
    {
        return -1;
    }

    // Find the nearest cell centre to this location
    label celli = findNearestCell(location);

    // If point is in the nearest cell return
    if (pointInCell(location, celli))
    {
        return celli;
    }
    else // point is not in the nearest cell so search all cells
    {
        bool cellFound = false;
        label n = 0;

        while ((!cellFound) && (n < nCells()))
        {
            if (pointInCell(location, n))
            {
                cellFound = true;
                celli = n;
            }
            else
            {
                n++;
            }
        }
        if (cellFound)
        {
            return celli;
        }
        else
        {
            return -1;
        }
    }
}


// ************************************************************************* //
