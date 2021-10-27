
#ifndef interTwoPhaseCentralFoam_H
#define interTwoPhaseCentralFoam_H

#include "fvm.H"
#include "fvc.H"
#include "uniformDimensionedFields.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarField.H"
#include "kappaFunction.H"
#include "standardMachToAcCourantRatioKappaFunction.H"
#include "psiThermo.H"
//#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class interTwoPhaseCentralFoam
:
   public IOdictionary
{
protected:
    //- Primitative variables

        //- Velocity
        volVectorField U_;

        //- Velocity estimated on cell surfaces
       surfaceScalarField phi_;
       surfaceScalarField phiC_;

        //- Pressure
        volScalarField p_;

        //- Temperature
        volScalarField T_;

        //Universal gas constant
        dimensioned< scalar > R_;

        //R/Molar mass
        dimensioned< scalar > R1_;
        dimensioned< scalar > R2_;

        //Molecular mass of the phase 1
        dimensioned< scalar > molM1_;

        //Molecular mass of the phase 2
        dimensioned< scalar > molM2_;

        //A coefficient that relates pressure and density of a phase 1
        volScalarField psi1_;

        //A coefficient that relates pressure and density of a phase 2
        volScalarField psi2_;

        //- Phase 1 density
        volScalarField rho1_;

        //- Phase 2 density
        volScalarField rho2_;


        //- Volume fraction of phase 1
        volScalarField volumeFraction1_;

        //- Volume fraction of phase 2
        volScalarField volumeFraction2_;

        //- Total mass
        volScalarField rho_;


    //- Constants

       //Constant density of the liquid
       dimensioned< scalar > rho01Sc_;
       dimensioned< scalar > rho02Sc_;

       volScalarField rho01_;
       volScalarField rho02_;

       // Specific heat (Cp) of the phase 1
       dimensioned< scalar > Cp1_;

       // Specific heat (Cp) of the phase 2
       dimensioned< scalar > Cp2_;

       //Dynamic viscousity of the phase 1
       dimensioned< scalar > mu1_;

       //Dynamic viscousity of the phase 2
       dimensioned< scalar > mu2_;

       //Dynamic viscousity of a mixture
       volScalarField mu_;

       //Kinematic viscousity of a mixture
       volScalarField nu_;

       dimensioned< scalar > Pr_;

       //Thermal conductivity of the phase 1
       dimensioned< scalar > lambda1_;

       //Thermal conductivity of the phase 2
       dimensioned< scalar > lambda2_;

       volScalarField lambda_;

       //Thermal diffusivity
       dimensioned< scalar > alpha1_;

       dimensioned< scalar > alpha2_;

       volScalarField alpha_;

       //Heat capacity ratio
       dimensioned< scalar > gamma1_;

       dimensioned< scalar > gamma2_;

      //Squared speed of sound
      volScalarField C1_;

      volScalarField C2_;

      volScalarField C_;

      //Acoustic impendans
      volScalarField Z1_;

      volScalarField Z2_;

      //A compressibiliy factor
      volScalarField K_;


      //The matrix H/A for the pressure equation
      volVectorField HbyA_;

      //The matrix 1/A for the pressure equation
      volScalarField rbyA_;

      //The magnitude of Velocity
      // it calculates: magSqr(U_)
      volScalarField rho1AmpU_;

      volScalarField rho2AmpU_;


/***********************Tadmor-Kurganov Scheme*******************************/
/*******************************Region One**********************************/

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

      surfaceScalarField phi1_;

/*******************************Region Two**********************************/

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

      surfaceScalarField phi2_;

/***********************Tadmor-Kurganov Scheme*******************************/

/*************************Pressure Equation*********************************/

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

      surfaceScalarField p_own;
      surfaceScalarField p_nei;

      volVectorField gradp1;
      volVectorField gradp2;

/************************kappa Function*********************************/
      autoPtr<fv::kappaFunction> kappaFuncPtr;

      surfaceScalarField kappa;

      surfaceScalarField amaxSf;

      surfaceScalarField uMagSf;

//      surfaceVectorField uSf;


//      dimensionedScalar amaxSmall;

/************************Viscosity Function*********************************/

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

     surfaceScalarField rho01_own;
     surfaceScalarField rho01_nei;

     surfaceScalarField rho02_own;
     surfaceScalarField rho02_nei;

     surfaceScalarField phi01d_own;
     surfaceScalarField phi01d_nei;

     surfaceScalarField phi02d_own;
     surfaceScalarField phi02d_nei;

     dimensionedScalar scAlpha_;
     surfaceScalarField phic_;

     surfaceScalarField vF1_;
     volScalarField E_;
     volScalarField E1_;
     volScalarField E2_;
     volScalarField EF_;

     volScalarField Ep_;


public:

   TypeName("interTwoPhaseCentralFoam");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interTwoPhaseCentralFoam,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );

// * * * * * * * * * *  Constructor and Destructor Functions * * * * * * * * //
    // Constructor
    interTwoPhaseCentralFoam(const fvMesh& mesh);

    // Selector
       static autoPtr<interTwoPhaseCentralFoam> New(const fvMesh& mesh);

    //- Destructor
    virtual ~interTwoPhaseCentralFoam();

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //Member Functions
    //volScalarField saveOld();
    void saveOld();

    //Return functions
    volScalarField rho1()
    {
      return rho1_;
    }

    surfaceScalarField phiC()
    {
      phiC_ = fvc::flux(U_);
      return phiC_;
    }
/*
    surfaceScalarField C()
    {
      C_ = fvc::interpolate(sqrt(C1_));
      C_.setOriented(true);
      return C_;
    }
*/
    volScalarField rho2()
    {
      return rho2_;
    }

    volScalarField rho()
    {
      return rho_;
    }

    volVectorField U()
    {
      return U_;
    }

    volScalarField T()
    {
      return T_;
    }

    volScalarField p()
    {
      return p_;
    }

    volScalarField volumeFraction1()
    {
      return volumeFraction1_;
    }

    volScalarField volumeFraction2()
    {
      return volumeFraction2_;
    }

    volScalarField mu()
    {
      return mu_;
    }

    //- Solve sub-step stepi
    virtual void alpha1Eqnsolve();
    virtual void alpha2Eqnsolve();
    virtual void UEqnsolve();
    virtual void TEqnsolve();
    virtual void rho1Eqnsolve();
    virtual void rho2Eqnsolve();
    virtual void pEqnsolve();

    //- Update
    virtual void Update();
    virtual void UpdateVelocity();
    virtual void Flux();
    virtual void alphaPhi();

//    surfaceScalarField rhoU1TKS_();
//    surfaceScalarField rhoU2TKS_();

    void rhoU1TKS_();
    void rhoU2TKS_();

    void Velocity();
    void DensityThermo();
    void DensityNull();
    void DensityCorrection();
    void Density();
    void ThermoCoefficient();
    void Initialize();
    void pressureGradient();
    void CallFlux();
    void UEqn();
    void UpdateCentralWeights();
    void UpdateCentralFields();

    void UpdateKappa();
    void divDevRhoReff();
    void Tviscosity();
    void alpha();
    void devRhoReff();
    void TViscousitySource();
    void info();
    void alphaPhiCorrect();


}; // End of the class


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
