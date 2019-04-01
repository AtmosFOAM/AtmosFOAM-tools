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
    
    // Check that levelInterfaces are monotonically increasing
    for(label il = 0; il < levelInterfaces.size()-1; il++)
    {
        if (levelInterfaces[il+1] <= levelInterfaces[il])
        {
            FatalErrorIn("horizontalMean")
                << "levelInterfaces must be monotonically increasing, not "
                << levelInterfaces << exit(FatalError);
        }
    }
    
    // Mid levels for output
    const label nLevels = levelInterfaces.size()-1;
    scalarList levels(nLevels);
    for(label il = 0; il < nLevels; il++)
    {
        levels[il] = 0.5*(levelInterfaces[il] + levelInterfaces[il+1]);
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // Find the total volume in all cells in each layer
        scalarList totalVolume(nLevels, scalar(0));
        forAll(mesh.V(), cellI)
        {
            const scalar h = mesh.C()[cellI].z();
            label ilh = -1;
            for(label il = 0; il < nLevels && ilh == -1; il++)
            {
                if (h >= levelInterfaces[il] && h <= levelInterfaces[il+1])
                {
                    ilh = il;
                }
            }
            
            if (ilh != -1)
            {
                totalVolume[ilh] += mesh.V()[cellI];
            }
        }

        forAll(fieldNamesAndDensities, fieldI)
        {
            const Pair<word>& fieldDensity = fieldNamesAndDensities[fieldI];
        
            const volScalarField f
            (
                IOobject
                (
                    fieldDensity[0], runTime.timeName(), mesh,
                    IOobject::MUST_READ
                ),
                mesh
            );
            const volScalarField rho
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

            // initialise output file
            fileName outFile = args.optionFound("cellSet") ?
                args.rootPath() / args.caseName() / runTime.timeName() 
                    / "horizontalMean_"+args.optionRead<string>("cellSet")
                      +"_"+fieldDensity[0]+".dat"
              : args.rootPath() / args.caseName() / runTime.timeName() 
                    / "horizontalMean_"+fieldDensity[0]+".dat";
            Info << "Writing horizontal means to " << outFile << endl;
            OFstream os(outFile);
            os << "#level volFraction rho mean stdDev min max" << endl;
            
            scalarList meanfRho(nLevels, scalar(0));
            scalarList meanSqrfRho(nLevels, scalar(0));
            scalarList minf(nLevels, GREAT);
            scalarList maxf(nLevels, -GREAT);
            scalarList mass(nLevels, scalar(0));
            scalarList volume(nLevels, scalar(0));
            
            // Loop through all cells and asign values to the correct level
            forAll(sumCells, i)
            {
                const label cellI = sumCells[i];
                const scalar h = mesh.C()[cellI].z();
                label ilh = -1;
                for(label il = 0; il < nLevels && ilh == -1; il++)
                {
                    if (h >= levelInterfaces[il] && h <= levelInterfaces[il+1])
                    {
                        ilh = il;
                    }
                }
                
                if (ilh != -1)
                {
                    volume[ilh] += mesh.V()[cellI];
                    mass[ilh] += mesh.V()[cellI]*rho[cellI];
                    meanfRho[ilh] += mesh.V()[cellI]*f[cellI]*rho[cellI];
                    meanSqrfRho[ilh] += mesh.V()[cellI]
                                       *sqr(f[cellI] )*rho[cellI];
                    if (f[cellI] <= minf[ilh]) minf[ilh] = f[cellI];
                    if (f[cellI] >= maxf[ilh]) maxf[ilh] = f[cellI];
                }
            }
            
            for(label il = 0; il < nLevels; il++)
            {
                scalar massMin = max(mass[il], VSMALL);
                scalar var = (massMin*meanSqrfRho[il] - sqr(meanfRho[il]))
                             /sqr(massMin);
                
                os << levels[il] << " " 
                   << volume[il]/totalVolume[il] << " " 
                   << massMin/volume[il] << " "
                   << meanfRho[il]/massMin << " "
                   << Foam::sqrt(var) << " " 
                   << minf[il] << " "
                   << maxf[il] << endl;
            }
        }
    }
    return(0);
}


// ************************************************************************* //
