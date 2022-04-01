/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    mapFields

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields time directories specified.
    Only for non-parallel cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToMesh0.H"
#include "processorFvPatch.H"
#include "MapMeshes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapConsistentMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
    else
    {
        MapConsistentMesh<eqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapSubMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            patchMap,
            cuttingPatches,
            mapOrder
        );
    }
    else
    {
        MapSubMesh<eqOp>
        (
            meshSource,
            meshTarget,
            patchMap,
            cuttingPatches,
            mapOrder
        );
    }
}


void mapConsistentSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapConsistentSubMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
    else
    {
        MapConsistentSubMesh<eqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
}


wordList addProcessorPatches
(
    const fvMesh& meshTarget,
    const wordList& cuttingPatches
)
{
    // Add the processor patches to the cutting list
    HashTable<label> cuttingPatchTable;
    forAll(cuttingPatches, i)
    {
        cuttingPatchTable.insert(cuttingPatches[i], i);
    }

    forAll(meshTarget.boundary(), patchi)
    {
        if (isA<processorFvPatch>(meshTarget.boundary()[patchi]))
        {
            if
            (
               !cuttingPatchTable.found
                (
                    meshTarget.boundaryMesh()[patchi].name()
                )
            )
            {
                cuttingPatchTable.insert
                (
                    meshTarget.boundaryMesh()[patchi].name(),
                    -1
                );
            }
        }
    }

    return cuttingPatchTable.toc();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::validArgs.append("sourceCase");

    timeSelector::addOptions();

    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
    argList::addBoolOption
    (
        "consistent",
        "source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "specify the mapping method"
    );
    argList::addBoolOption
    (
        "subtract",
        "subtract mapped source from target"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = fvMesh::defaultRegion;
    if (args.optionFound("sourceRegion"))
    {
        sourceRegion = args["sourceRegion"];
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;
    if (args.optionFound("targetRegion"))
    {
        targetRegion = args["targetRegion"];
        Info<< "Target region: " << targetRegion << endl;
    }

    const bool consistent = args.optionFound("consistent");

    meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;
    if (args.optionFound("mapMethod"))
    {
        const word mapMethod(args["mapMethod"]);
        if (mapMethod == "mapNearest")
        {
            mapOrder = meshToMesh0::MAP;
        }
        else if (mapMethod == "interpolate")
        {
            mapOrder = meshToMesh0::INTERPOLATE;
        }
        else if (mapMethod == "cellPointInterpolate")
        {
            mapOrder = meshToMesh0::CELL_POINT_INTERPOLATE;
        }
        else if (mapMethod == "cellVolumeWeight")
        {
            mapOrder = meshToMesh0::CELL_VOLUME_WEIGHT;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mapMethod " << mapMethod << ". Valid options are: "
                << "cellVolumeWeight, mapNearest, interpolate and cellPointInterpolate"
                << exit(FatalError);
        }

        Info<< "Mapping method: " << mapMethod << endl;
    }

    const bool subtract = args.optionFound("subtract");
    if (subtract)
    {
        Info<< "Subtracting mapped source field from target" << endl;
    }

    #include "createTimes.H"
    instantList timeDirs = timeSelector::select0(runTimeSource, args);

    HashTable<word> patchMap;
    wordList cuttingPatches;

    if (!consistent)
    {
        IOdictionary mapFieldsDict
        (
            IOobject
            (
                "mapFieldsDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        mapFieldsDict.lookup("patchMap") >> patchMap;
        mapFieldsDict.lookup("cuttingPatches") >>  cuttingPatches;
    }

    // Mapping for serial source and target for all times

    Info<< "Create meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTimeSource.setTime(timeDirs[timeI], timeI);
        runTimeTarget.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTimeSource.timeName() << endl;

        meshSource.readUpdate();
        meshTarget.readUpdate();

        if (consistent)
        {
            mapConsistentMesh(meshSource, meshTarget, mapOrder, subtract);
        }
        else
        {
            mapSubMesh
            (
                meshSource,
                meshTarget,
                patchMap,
                cuttingPatches,
                mapOrder,
                subtract
            );
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
