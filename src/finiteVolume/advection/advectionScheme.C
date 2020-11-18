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

\*---------------------------------------------------------------------------*/

#include "advectionScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatchField.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::advectionScheme<Type>>
Foam::advectionScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshFluxConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (debug)
    {
        InfoInFunction
            << "Discretisation scheme = " << schemeName << endl;
    }

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshFluxConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshFluxConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::advectionScheme<Type>::~advectionScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::advectionScheme<Type>::advect
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& faceFlux
) const
{
    if (debug)
    {
        InfoInFunction
            << "Advecting "
            << vf.type() << " "
            << vf.name() << " with fluxes " << faceFlux.name()
            << endl;
    }

    upwind<Type> uInterp(vf.mesh(), faceFlux);
    
    tmp<GeometricField<Type, fvPatchField, volMesh>> tdiv
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            "advetion("+faceFlux.name()+vf.name()+')',
            fvc::surfaceIntegrate(faceFlux*uInterp.interpolate(vf))
        )
    );

    return tdiv;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::advectionScheme<Type>::advect
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const surfaceScalarField& faceFlux
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tdiv
        = advect(tvf(), faceFlux);
    tvf.clear();
    return tdiv;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::advectionScheme<Type>::advect
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tfaceFlux
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tdiv
        = advect(vf, tfaceFlux());
    tfaceFlux.clear();
    return tdiv;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::advectionScheme<Type>::advect
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const tmp<surfaceScalarField>& tfaceFlux
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tdiv
        = advect(tvf(), tfaceFlux());
    tvf.clear();
    tfaceFlux.clear();
    return tdiv;
}


// ************************************************************************* //
