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

Class
    Foam::mulesWithDiffusionImplicitLimiter

Description
    MULES limiter adapted for using with hybrid Kurganov-Tadmor scheme

SourceFiles
    customMULES.C

\*---------------------------------------------------------------------------*/

#ifndef customMULES_H
#define customMULES_H

#include "gaussConvectionScheme.H"
#include "MULES.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "scalarField.H"
#include "slicedSurfaceFieldsFwd.H"
#include "PtrList.H"


namespace Foam
{
    void mulesWithDiffusionImplicitLimiter
    (
        const volScalarField& rho,
        volScalarField& Y,
        const surfaceScalarField& phi_own,
        const surfaceScalarField& phi_nei,
        scalarField& lambdaCoeff,
        surfaceScalarField& rhoPhifCorr, //total flux as input, correction flux as output
        surfaceScalarField& diffFlux,
        const surfaceScalarField& Dm,
        const fvScalarMatrix& SuSp
    );

    void mulesWithDiffusionImplicitLimiter
    (
        const volScalarField& rho,
        volScalarField& Y,
        const surfaceScalarField& phi,
        scalarField& lambdaCoeff,
        surfaceScalarField& rhoPhifCorr, //total flux as input, correction flux as output
        surfaceScalarField& diffFlux,
        const surfaceScalarField& Dm,
        const fvScalarMatrix& SuSp
    );
};

#endif
