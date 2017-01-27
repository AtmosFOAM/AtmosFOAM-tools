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
    setGaussians.C

Description
    Initial a given volScalarField to be the sum of a set of Gaussians with 
    different maxima, radius and centres

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gaussian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("dictionary name (in system)");
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"
#   include "createTime.H"
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
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word dictName = args.args()[1].c_str();
    Info<< "Reading initial conditions from" << dictName << endl;

    IOdictionary initDict
    (
        IOobject
        (
            dictName, mesh.time().system(), mesh, IOobject::MUST_READ
        )
    );
    
    // Read in background value
    dimensionedScalar value
    (
        initDict.lookupOrDefault<dimensionedScalar>
        (
            "backgroundValue",
            dimensionedScalar("T", dimless, scalar(0))
        )
    );

    // Read in list of Gaussians
    List<gaussian> gaussians(initDict.lookup("gaussians"));

    // Initialise tracer to the background value
    volScalarField T
    (
        IOobject(value.name(), runTime.timeName(), mesh),
        mesh,
        value
    );

    // Check that all Gaussians and the background value have the same
    // dimensions
    forAll(gaussians, ig)
    {
        if (gaussians[ig].max().dimensions() != value.dimensions())
        {
            FatalErrorIn("setGaussians")
                << " backgroundValue defined with dimensions "
                << value.dimensions() << " but Gaussian[" << ig
                << "] max defined with dimensions "
                << gaussians[ig].max().dimensions() << exit(FatalError);
        }
    }
    
    // Add fields for all of the Gaussian distributions
    forAll(gaussians, ig)
    {
        T += gaussians[ig].field(mesh);
    }
    
    T.write();
    
    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
