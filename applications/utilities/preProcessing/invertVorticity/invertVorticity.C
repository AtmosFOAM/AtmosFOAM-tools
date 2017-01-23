/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Application
    invertVorticity.C

Description
    Reads in the scalar valued vorticity on a 2d mesh and inverts to find the
    streamfunciton and velocity. Needs to read in the velocity field, Uf for
    the boundary conditions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvcCurlf.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("dictionary name (in system)");
    timeSelector::addOptions();
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word dictName = args.args()[1].c_str();
    Info << "Read in constant background wind from" << dictName << endl;
    IOdictionary initDict
    (
        IOobject
        (
            dictName, mesh.time().system(), mesh, IOobject::MUST_READ
        )
    );
    const dimensionedVector U0(initDict.lookup("U0"));

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        Info << "Mesh has normal direction" << flush;
        const vector meshNormal = 0.5*(Vector<label>(1,1,1)-mesh.geometricD());
        Info << meshNormal << endl;

        Info << "Reading in vorticity" << endl;
        volScalarField vorticity
        (
            IOobject("vorticity", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );

        Info << "Reading in velocity field on faces, Uf" << endl;
        surfaceVectorField Uf
        (
            IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        
        Info << "Creating the streamFunction" << endl;
        volScalarField streamFunction
        (
            IOobject("streamFunction", runTime.timeName(), mesh,
                     IOobject::MUST_READ),
            mesh
        );
        
        // Set the streamfunction for the background uniform flow
        const dimensionedVector velocityPerp = meshNormal ^ U0;
        streamFunction == magSqr(U0)/magSqr(velocityPerp)
                         *(mesh.C() & velocityPerp);

        // Invert the streamfunction to find the vorticity
        bool converged = false;
        for(label it = 0; it < 100 && !converged; it++)
        {
            fvScalarMatrix streamFuncEqn
            (
                fvm::laplacian(streamFunction) == vorticity
            );
            solverPerformance sp = streamFuncEqn.solve();
            converged = sp.nIterations() <= 0;
        }

        // Set the velocity from the streamfunction
        Uf = linearInterpolate(fvc::curl(streamFunction*meshNormal));

        Uf.write();
        streamFunction.write();
    }
    
    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
