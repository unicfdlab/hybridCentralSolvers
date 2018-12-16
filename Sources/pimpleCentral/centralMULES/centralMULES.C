/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
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

#include "centralMULES.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "slicedSurfaceFields.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvm.H"
#include "gaussConvectionScheme.H"
#include "fvMatrices.H"


void Foam::mulesWithDiffusionImplicitLimiter
(
    const volScalarField& rho,
    volScalarField& Y,
    const surfaceScalarField& phi_own,
    const surfaceScalarField& phi_nei,
    scalarField& lambdaFace,
    surfaceScalarField& rhoPhif,
    surfaceScalarField& diffFlux,
    const surfaceScalarField& Dmi,
    const fvScalarMatrix& SuSp
)
{
    const fvMesh& mesh = rho.mesh();
    const word Yname (Y.name());
    Y.rename("Yi");

    upwind<scalar> UDsOwn(mesh, phi_own);
    upwind<scalar> UDsNei(mesh, phi_nei);

    fvScalarMatrix YConvection
    (
        fv::gaussConvectionScheme<scalar>(mesh, phi_own, UDsOwn).fvmDiv(phi_own, Y)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_nei, UDsNei).fvmDiv(phi_nei, Y)
    );

    surfaceScalarField rhoPhifBD = YConvection.flux();

    surfaceScalarField& rhoPhifCorr = rhoPhif;
    rhoPhifCorr -= rhoPhifBD;

    volScalarField Su
    (
        "Su",
        SuSp & Y
    );

    MULES::limiter
    (
        lambdaFace,
        1.0/mesh.time().deltaTValue(),
        rho,
        Y,
        rhoPhifBD,
        rhoPhifCorr,
        zeroField(), //Sp
        Su,
        1.0, //psiMax,
        0.0 //psiMin,
    );

    Y.rename(Yname);
}

void Foam::mulesWithDiffusionImplicitLimiter
(
    const volScalarField& rho,
    volScalarField& Y,
    const surfaceScalarField& phi,
    scalarField& lambdaFace,
    surfaceScalarField& rhoPhif,
    surfaceScalarField& diffFlux,
    const surfaceScalarField& Dmi,
    const fvScalarMatrix& SuSp
)
{
    const fvMesh& mesh = rho.mesh();
    const word Yname (Y.name());
    Y.rename("Yi");

    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix YConvection
    (
        fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, Y)
    );

    surfaceScalarField rhoPhifBD = YConvection.flux();

    surfaceScalarField& rhoPhifCorr = rhoPhif;
    rhoPhifCorr -= rhoPhifBD;

    volScalarField Su
    (
        "Su",
        SuSp & Y
    );

    MULES::limiter
    (
        lambdaFace,
        1.0/mesh.time().deltaTValue(),
        rho,
        Y,
        rhoPhifBD,
        rhoPhifCorr,
        zeroField(),
        Su,
        1.0, //psiMax,
        0.0 //psiMin,
    );
    Y.rename(Yname);
}

//
//END-OF-FILE
//
