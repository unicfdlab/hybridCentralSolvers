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
	const surfaceScalarField& phiPos,
	const surfaceScalarField& phiNeg,
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

