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
    const surfaceScalarField& phiPos,
    const surfaceScalarField& phiNeg,
    scalarField& lambdaFace,
    surfaceScalarField& rhoPhif,
    surfaceScalarField& diffFlux,
    const surfaceScalarField& Dmi,
    const fvScalarMatrix& SuSp
)
{
    const fvMesh& mesh = rho.mesh();
    const dictionary& MULEScontrols = mesh.solutionDict().solverDict("Yi");

    label nLimiterIter
    (
	readLabel(MULEScontrols.lookup("nLimiterIter"))
    );
    
    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);
    
    fvScalarMatrix YConvection
    (
	fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).fvmDiv(phiPos, Y)
	+
	fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsPos).fvmDiv(phiNeg, Y)
    );
    
    surfaceScalarField rhoPhifBD = YConvection.flux();

    surfaceScalarField& rhoPhifCorr = rhoPhif;
    rhoPhifCorr -= rhoPhifBD;

    volScalarField Su
    (
	"Su",
	SuSp & Y
    );

//1.0/mesh.time().deltaTValue(),

    MULES::limiter
    (
	lambdaFace,
	rho,
	Y,
	rhoPhifBD,
	rhoPhifCorr,
	zeroField(),
	Su,
	1.0, //psiMax,
	0.0,  //psiMin,
	nLimiterIter
    );
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
    const dictionary& MULEScontrols = mesh.solutionDict().solverDict("Yi");

    label nLimiterIter
    (
	readLabel(MULEScontrols.lookup("nLimiterIter"))
    );
    
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
	rho,
	Y,
	rhoPhifBD,
	rhoPhifCorr,
	zeroField(),
	Su,
	1.0, //psiMax,
	0.0, //psiMin,
	nLimiterIter
    );
}

//
//END-OF-FILE
//

