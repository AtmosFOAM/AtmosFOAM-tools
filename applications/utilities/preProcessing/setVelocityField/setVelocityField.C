#include "fvCFD.H"
#include "velocityField.H"

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "dict", "dictName", "specify the dictionary name (in system)"
    );
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
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

    Info << "Creating flux field phi" << endl;
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("phi", cmptMultiply(dimVelocity, dimArea), scalar(0)),
       "fixedValue"
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::reconstruct(phi)
    );

    // Read Uf if present, otherwise create and write (not used)
    surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U)
    );

    const word dictName = args.optionFound("dict") ?
                          args.optionRead<word>("dict") :
                          "velocityFieldDict";
    Info<< "Reading initial conditions from " << dictName << endl;
    IOdictionary dict
    (
        IOobject
        (
            dictName,
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<velocityField> v(velocityField::New(dict));

    forAll(timeDirs, timeI)
    {
        Info << "writing phi for time " << runTime.timeName() << endl;

        v->applyTo(phi);
        phi.write();

        U = fvc::reconstruct(phi);
        Uf = linearInterpolate(U);
        Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        U.write();
        Uf.write();
    }

    return EXIT_SUCCESS;
}


