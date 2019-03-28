/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    horizontalMean

Description
    Calculates the horizontal mean, standard deviation, min and max for
    a list of variables given in dictionary system/horizontalMeanDict.
    Also in the dictionary is a list of vertical level boundaries and for each
    field, need to say the density variable used to calculate volume average.
    Can use a subset of cells from a cellSet

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "OFstream.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    argList::addOption
    (
        "cellSet", "cellSetName", "only calculate sums for a subset of cells"
    );

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // List of cells to sum
    labelList sumCells;
    if (args.optionFound("cellSet"))
    {
        const word cellSetName(args.optionRead<string>("cellSet"));
        cellSet cells(mesh, cellSetName);
        sumCells = cells.toc();
    }
    else
    {
        // Select all cells
        sumCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            sumCells[cellI] = cellI;
        }
    }

    Info << "\nReading dictionary horizontalMeanDict\n" << endl;
    IOdictionary hMeanDict
    (
        IOobject
        (
            "horizontalMeanDict", runTime.system(), runTime,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    
    const scalarList levelInterfaces(hMeanDict.lookup("levelInterfaces"));
    const List<Pair<word>> fieldNamesAndDensities
    (
        hMeanDict.lookup("fieldNamesAndDensities")
    );
    
    // Mid levels for output
    scalarList levels(levelInterfaces.size()-1);
    for(label il = 0; il < levels.size(); il++)
    {
        levels[il] = 0.5*(levelInterfaces[il] + levelInterfaces[il+1]);
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        forAll(fieldNamesAndDensities, fieldI)
        {
            const Pair<word>& fieldDensity = fieldNamesAndDensities[fieldI];
            Info << "Reading in field and density pair " << fieldDensity << endl;
            

            // initialise output file
            fileName outFile = args.optionFound("cellSet") ?
                args.rootPath() / args.caseName() / runTime.timeName() 
                    / "horizontalMean"+args.optionRead<string>("cellSet")+fieldDensity[0]+".dat"
              : args.rootPath() / args.caseName() / runTime.timeName() 
                    / "horizontalMean"+fieldDensity[0]+".dat";
            Info << "Writing horizontal means to " << outFile << endl;
            OFstream os(outFile);
            os << "#level mean stdDev min max" << endl;
        
            volScalarField f
            (
                IOobject
                (
                    fieldDensity[0], runTime.timeName(), mesh,
                    IOobject::MUST_READ
                ),
                mesh
            );
            volScalarField rho
            (
                IOobject
                (
                    fieldDensity[1], runTime.timeName(), mesh,
                    IOobject::READ_IF_PRESENT
                ),
                mesh,
                dimensionedScalar("one", dimless, scalar(1))
            );
            Info << fieldDensity[0] << " weighted by density " 
                 << fieldDensity[1] << endl;
        }
    }
    return(0);
}


// ************************************************************************* //
