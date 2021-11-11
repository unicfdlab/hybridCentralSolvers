
#ifndef interTwoPhaseCentralFoam_H
#define interTwoPhaseCentralFoam_H

#include "fvm.H"
#include "fvc.H"
#include "uniformDimensionedFields.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarField.H"
//#include "kappaFunction.H"
//#include "standardMachToAcCourantRatioKappaFunction.H"
#include "psiThermo.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interTwoPhaseCentralFoam Declaration
\*---------------------------------------------------------------------------*/

class interTwoPhaseCentralFoam
:
   public IOdictionary
{
protected:

        pimpleControl& pimple_;

        //- Velocity
        //- MUST_READ
        volVectorField U_;

        //- Velocity estimated on cell boundaries
        surfaceScalarField phi_;

        //- Pressure
        //- MUST_READ
        volScalarField p_;

        //- Temperature
        //- MUST_READ
        volScalarField T_;

        //- Universal gas constant
        //- Set in the phaseProperties file
        dimensioned< scalar > R_;

        //- Specific gas constant of the first phase
        //- Set in the phaseProperties file
        dimensioned< scalar > R1_;

        //- Specific gas constant of the second phase
        //- Set in the phaseProperties file
        dimensioned< scalar > R2_;

        //- Molecular mass of the first phase
        //- Set in the phaseProperties file
        dimensioned< scalar > molM1_;

        //- Molecular mass of the second phase
        //- Set in the phaseProperties file
        dimensioned< scalar > molM2_;

        //- Compressibility coefficient of the first phase
        volScalarField psi1_;

        //- Compressibility coefficient of the second phase
        volScalarField psi2_;

        //- Density of the first phase
        volScalarField rho1_;

        //- Density of the second phase
        volScalarField rho2_;


        //- Volume fraction of the first phase
        //- MUST_READ
        volScalarField volumeFraction1_;

        //- Volume fraction of the first phase
        //- NO_READ
        volScalarField volumeFraction2_;

        //- Total density of a mixture
        volScalarField rho_;

       // Constant density of the first phase, i.e. rho_01
       // Set in the phaseProperties file as rho01
       dimensioned< scalar > rho01Sc_;

       // Constant density of a phase, i.e. rho_02
       // Set in the phaseProperties file as rho02
       dimensioned< scalar > rho02Sc_;

       //- Constant density of the first phase, i.e. rho_01
      //- Set in the phaseProperties file as rho01
       volScalarField rho01_;

       //- Constant density of the first phase, i.e. rho_02
       //- Set in the phaseProperties file as rho02
       volScalarField rho02_;

       //- Specific heat of the first phase
       //- Set in the phaseProperties
       dimensioned< scalar > Cp1_;

       //- Specific heat of the second phase
       //- Set in the phaseProperties
       dimensioned< scalar > Cp2_;

       //- Dynamic viscousity of the first phase
       //- Set in the phaseProperties
       dimensioned< scalar > mu1_;

       //- Dynamic viscousity of the second phase
       //- Set in the phaseProperties
       dimensioned< scalar > mu2_;

       //- Prandtl number
       //- Set in the phaseProperties
       dimensioned< scalar > Pr_;

       //- Thermal diffusivity of the first phase
       dimensioned< scalar > alpha1_;

       //- Thermal diffusivity of the second phase
       dimensioned< scalar > alpha2_;

       //- Heat capacity ratio of the first phase
       dimensioned< scalar > gamma1_;

       //- Heat capacity ratio of the second phase
       dimensioned< scalar > gamma2_;

      //- Squared speed of sound of the first phase
      volScalarField C1_;

      //- Squared speed of sound of the second phase
      volScalarField C2_;

      //- Squared speed of sound of a mixture
      volScalarField C_;

      //- Acoustic impendans of the first phase
      volScalarField Z1_;

      //- Acoustic impendans of the second phase
      volScalarField Z2_;

      //- A compressibiliy factor
      volScalarField K_;

      //- The H/A matrix of the Momentum equation
      volVectorField HbyA_;

      //- The 1/A matrix of the Momentum equation
      volScalarField rbyA_;

      //- The magnitude of Velocity of the first phase
      // It is calculated as rho1_*magSqr(U_)
      volScalarField rho1AmpU_;

      //- The magnitude of Velocity of the second phase
      // It is calculated as rho2_*magSqr(U_)
      volScalarField rho2AmpU_;

/********************Variables for Update-Central Wieghts Function*************/

      dimensionedScalar v_zero;

      surfaceScalarField own_;

      surfaceScalarField nei_;

      surfaceScalarField rho1_own;

      surfaceScalarField rho1_nei;

      surfaceScalarField phi1v_own;

      surfaceScalarField phi1v_nei;

      surfaceScalarField C1_own;

      surfaceScalarField C1_nei;

      surfaceScalarField C1Sf_own;

      surfaceScalarField C1Sf_nei;

      surfaceScalarField ap1;

      surfaceScalarField am1;

      surfaceScalarField alpha1_own;

      surfaceScalarField aSf1;

      surfaceScalarField alpha1_nei;

      surfaceScalarField phi1_own;

      surfaceScalarField phi1_nei;

/**********************************Phase Two***********************************/

      surfaceScalarField rho2_own;

      surfaceScalarField rho2_nei;

      surfaceScalarField phi2v_own;

      surfaceScalarField phi2v_nei;

      surfaceScalarField C2_own;

      surfaceScalarField C2_nei;

      surfaceScalarField C2Sf_own;

      surfaceScalarField C2Sf_nei;

      surfaceScalarField ap2;

      surfaceScalarField am2;

      surfaceScalarField alpha2_own;

      surfaceScalarField aSf2;

      surfaceScalarField alpha2_nei;

      surfaceScalarField phi2_own;

      surfaceScalarField phi2_nei;

/****************************Update Central Fields*****************************/

      surfaceScalarField psi1_own;

      surfaceScalarField psi1_nei;

      surfaceVectorField psiU1_own;

      surfaceVectorField psiU1_nei;

      surfaceScalarField aphi1v_own;

      surfaceScalarField aphi1v_nei;

      surfaceScalarField phi1d_own;

      surfaceScalarField phi1d_nei;

      surfaceScalarField Dp1_own;

      surfaceScalarField Dp1_nei;

      surfaceScalarField psi2_own;

      surfaceScalarField psi2_nei;

      surfaceVectorField psiU2_own;

      surfaceVectorField psiU2_nei;

      surfaceScalarField aphi2v_own;

      surfaceScalarField aphi2v_nei;

      surfaceScalarField phi2d_own;

      surfaceScalarField phi2d_nei;

      surfaceScalarField Dp2_own;

      surfaceScalarField Dp2_nei;

      fvScalarMatrix pEqn1_own;

      fvScalarMatrix pEqn1_nei;

      fvScalarMatrix pEqn2_own;

      fvScalarMatrix pEqn2_nei;

/*******************************Reconstruct Velocity***************************/

      surfaceScalarField p_own;

      surfaceScalarField p_nei;

      volVectorField gradp1;

      volVectorField gradp2;

/******************************Viscosity Terms*********************************/

     fvVectorMatrix divDevRhoReff1_;

     fvVectorMatrix divDevRhoReff2_;

     fvScalarMatrix Tviscosity1;

     fvScalarMatrix Tviscosity2;

     volSymmTensorField devRhoReff1_;

     volSymmTensorField devRhoReff2_;

     volScalarField TSource1_;

     volScalarField TSource2_;

     surfaceScalarField psi1f;

     surfaceScalarField psi2f;

     surfaceScalarField phi01d_own;

     surfaceScalarField phi01d_nei;

     surfaceScalarField phi02d_own;

     surfaceScalarField phi02d_nei;

/******************************Something Essential*****************************/

     //- Continuity error of the first phase
     volScalarField E1_;

     //- Continuity error of the first phase
     volScalarField E2_;

     //- Continuity error of a mixture
     volScalarField E_;

     volScalarField ddtvF1_;

     volScalarField ddtvF2_;

     surfaceScalarField vF1face_;

     surfaceScalarField vF2face_;

     surfaceScalarField vFPhi1_;

     surfaceScalarField vFPhi2_;

     surfaceScalarField phiU_own;

     surfaceScalarField phiU_nei;

public:

   TypeName("interTwoPhaseCentralFoam");

// * * * * * * * * * *  Constructor and Destructor Functions * * * * * * * * //
    // Constructor
    interTwoPhaseCentralFoam(const fvMesh& mesh, pimpleControl& ctrl);

    //- Destructor
    virtual ~interTwoPhaseCentralFoam();

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void saveOld();


    surfaceScalarField phiC()
    {
      surfaceScalarField phiC = fvc::flux(U_);
      return phiC;
    }

// * * * * * * * * * * * * * * * * Main Functions  * * * * * * * * * * * * * //

    //- Solve an equation for volumeFraction1_
    virtual void alpha1Eqnsolve();

    //- Construct Momentum equation and extract rbyA_ and HbyA_
    void UEqn();

    //- Solve an equation for Temperature
    virtual void TEqnsolve();

    //- Estimate error of Continuity equation for the first phase
    virtual void rho1Eqnsolve();

    //- Estimate error of Continuity equation for the second phase
    virtual void rho2Eqnsolve();

    //- Solve an equation for pressure
    virtual void pEqnsolve();

// * * * * * * * * * * * * * * * * Middle Functions * * * * * * * * * * * * * //

    virtual void Flux();

    virtual void CompressibilityCoefficient();

    virtual void volumeFlux();

    void ReconstructVelocity();

    void DensityThermo();

    void Density();

    void Compressibility();

    void Initialize();

    void pressureGradient();

    void UpdateCentralWeights();

    void UpdateCentralFields();

    void divDevRhoReff();

    void Tviscosity();

    void alpha();

    void devRhoReff();

    void TViscousitySource();

    void info();

    void divU();

    void speedOfSound();

}; // End of the class


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
