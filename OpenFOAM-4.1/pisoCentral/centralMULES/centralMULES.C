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
        0.0  //psiMin,
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
        0.0  //psiMin,
    );
    Y.rename(Yname);
}

//
//END-OF-FILE
//

