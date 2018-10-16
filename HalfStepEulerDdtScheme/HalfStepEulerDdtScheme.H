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


    \*---------------------------------------------------------------------------*/

#ifndef HalfStepEulerDdtScheme_H
#define HalfStepEulerDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
    Class HalfStepEulerDdtScheme Declaration
    \*---------------------------------------------------------------------------*/

template<class Type>
class HalfStepEulerDdtScheme
:
    public ddtScheme<Type>
{
    // Private Member Functions

    //- No copy construct
    HalfStepEulerDdtScheme(const HalfStepEulerDdtScheme&) = delete;

    //- No copy assignment
    void operator=(const HalfStepEulerDdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("HalfStepEuler");


    // Constructors

    //- Construct from mesh
    HalfStepEulerDdtScheme(const fvMesh& mesh)
        :
        ddtScheme<Type>(mesh)
        {}

    //- Construct from mesh and Istream
    HalfStepEulerDdtScheme(const fvMesh& mesh, Istream& is)
        :
        ddtScheme<Type>(mesh, is)
        {}


    // Member Functions

    //- Return mesh reference
    const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

    tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
    (
        const dimensioned<Type>&
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
    (
        const dimensionedScalar&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
    (
        const volScalarField&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& psi
    );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fvcDdt
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    tmp<fvMatrix<Type>> fvmDdt
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    
    tmp<fvMatrix<Type>> fvmDdt
    (
        const dimensionedScalar&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    tmp<fvMatrix<Type>> fvmDdt
    (
        const volScalarField&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    tmp<fvMatrix<Type>> fvmDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& psi
    );

    typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

    tmp<fluxFieldType> fvcDdtUfCorr
    (
        const GeometricField<Type, fvPatchField, volMesh>& U,
        const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
    );

    tmp<fluxFieldType> fvcDdtPhiCorr
    (
        const GeometricField<Type, fvPatchField, volMesh>& U,
        const fluxFieldType& phi
    );

    tmp<fluxFieldType> fvcDdtUfCorr
    (
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& U,
        const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
    );

    tmp<fluxFieldType> fvcDdtPhiCorr
    (
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& U,
        const fluxFieldType& phi
    );

    tmp<surfaceScalarField> meshPhi
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );
};


template<>
tmp<surfaceScalarField> HalfStepEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> HalfStepEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> HalfStepEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> HalfStepEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "HalfStepEulerDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
