#include "interTwoPhaseCentralFoam.H"

namespace Foam
{
    defineTypeNameAndDebug(interTwoPhaseCentralFoam, 0);
    defineRunTimeSelectionTable(interTwoPhaseCentralFoam, dictionary);
}

Foam::interTwoPhaseCentralFoam::interTwoPhaseCentralFoam(const fvMesh& mesh, pimpleControl& ctrl)
:

    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    pimple_ (ctrl),

    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi_
    (
        "phi",
        (fvc::interpolate(U_))&mesh.Sf()
    ),

    phiC_
    (
        (fvc::interpolate(U_))&mesh.Sf()
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    R_
    (
      dimensioned< scalar >("R", *this)
    ),

    R1_
    (
      dimensioned< scalar >("R1", *this)
    ),

    R2_
    (
      dimensioned< scalar >("R2", *this)
    ),

    molM1_
    (
      dimensioned< scalar >("molM1", *this)
    ),

    molM2_
    (
      dimensioned< scalar >("molM2", *this)
    ),

    psi1_
    (
      "psi1",
      1 / (R1_ * T_)
    ),

    psi2_
    (
      "psi2",
      1 / (R2_ * T_)
    ),

    rho1_
    (
      IOobject
      (
          "rho1",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      psi1_*p_
    ),

    rho2_
    (
      IOobject
      (
          "rho2",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      psi2_*p_
    ),

    volumeFraction1_
    (
        IOobject
        (
            "volumeFraction1",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    volumeFraction2_
    (
      "volumeFraction2",
      1 - volumeFraction1_
    ),
/*
    volumeFraction2_
    (
      IOobject
      (
          "volumeFraction2",
          mesh.time().timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
    ),
*/
    rho_
    (
      IOobject
      (
          "rho",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      volumeFraction1_*rho1_ + volumeFraction2_*rho2_
    ),

    rho01Sc_
    (
      dimensioned< scalar >("rho01", *this)
    ),

    rho02Sc_
    (
      dimensioned< scalar >("rho02", *this)
    ),

    rho01_
    (
      IOobject
      (
          "rho2",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      (rho01Sc_*rho1_)/rho1_
    ),

    rho02_
    (
      IOobject
      (
          "rho2",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      (rho02Sc_*rho2_)/rho2_
    ),

    Cp1_
    (
      dimensioned< scalar >("Cp1", *this)
    ),

    Cp2_
    (
      dimensioned< scalar >("Cp2", *this)
    ),

    mu1_
    (
      dimensioned< scalar >("mu1", *this)
    ),

    mu2_
    (
      dimensioned< scalar >("mu2", *this)
    ),

    mu_
    (
      "mu",
      volumeFraction1_*mu1_ + volumeFraction2_*mu2_
      /*
      ((volumeFraction1_ * rho1_ + volumeFraction2_ * rho2_ )
      /
      (volumeFraction1_ * (rho1_ / mu1_) + volumeFraction2_ * (rho2_ / mu2_)))
      */
    ),

    nu_
    (
      "nu",
      mu_/rho_
    ),

    Pr_
    (
      dimensioned< scalar >("Pr", *this)
    ),

    lambda1_
    (
      "lambda1",
      (Cp2_*mu1_)/Pr_
    ),

    lambda2_
    (
      "lambda2",
      (Cp2_*mu1_)/Pr_
    ),

    lambda_
    (
      "lambda",
      volumeFraction1_*lambda1_ + volumeFraction2_*lambda2_
    ),

    alpha1_
    (
      "alpha1",
      mu1_/(Pr_)
    ),

    alpha2_
    (
      "alpha2",
      mu2_/(Pr_)
    ),

    alpha_
    (
      "alpha",
      volumeFraction1_*alpha1_+volumeFraction2_*alpha2_
    ),

    gamma1_
    (
      dimensioned< scalar >("gamma1", Cp1_/(Cp1_-(R_/molM1_)))
    ),

    gamma2_
    (
      dimensioned< scalar >("gamma2", Cp2_/(Cp2_-(R_/molM2_)))
    ),

    C1_
    (
      "C1",
      gamma1_*R1_*T_
    ),

    C2_
    (
      "C2",
      gamma2_*R2_*T_
    ),

    C_
    (
      "C",
      volumeFraction1_*sqrt(C1_) + volumeFraction2_*sqrt(C2_)
    ),

    Z1_
    (
      "Z1",
      rho1_*C1_
    ),

    Z2_
    (
      "Z2",
      rho2_*C2_
    ),

    K_
    (
      "K",
      (volumeFraction1_*volumeFraction2_*(Z1_-Z2_))
      /
      (Z1_*volumeFraction2_+Z2_*volumeFraction1_)
    ),

    HbyA_
    (
        IOobject
        (
            "HbyA",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 1, -1, 0, 0, 0, 0)
    ),

    rbyA_
    (
        IOobject
        (
            "rbyA",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(-1, 3, 1, 0, 0, 0, 0)
    ),

    rho1AmpU_
    (
      "rho1AmplitudeU",
      magSqr(U_) * rho1_
    ),

    rho2AmpU_
    (
      "rho2AmplitudeU",
      magSqr(U_) * rho2_
    ),


/***********************Tadmor-Kurganov Scheme*******************************/

    v_zero
    (
      "v_zero",
      dimVolume/dimTime,
      0.0
    ),

    own_
    (
        IOobject
        (
            "own",
            mesh.time().timeName(),
            //runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    ),

    nei_
    (
        IOobject
        (
            "nei",
            mesh.time().timeName(),
            //runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    ),

    rho1_own
    (
        "rho1own",
        fvc::interpolate(rho1_, own_, "reconstruct(rho1)")
    ),

    rho1_nei
    (
        "rho1nei",
        fvc::interpolate(rho1_, nei_, "reconstruct(rho1)")
    ),

    phi1v_own
    (
        "phi1v_own",
        (fvc::interpolate(rho1_*U_, own_, "reconstruct(U)") / rho1_own) & mesh.Sf()
    ),

    phi1v_nei
    (
        "phi1v_nei",
        (fvc::interpolate(rho1_*U_, nei_, "reconstruct(U)") / rho1_nei) & mesh.Sf()
    ),

     C1_own
    (
        "C1_own",
        fvc::interpolate(sqrt(C1_), own_, "reconstruct(phi)")
    ),

     C1_nei
    (
        "C1_nei",
        fvc::interpolate(sqrt(C1_), nei_, "reconstruct(phi)")
    ),

     C1Sf_own
    (
        "C1Sf_own",
        C1_own * mesh.magSf()
    ),

     C1Sf_nei
    (
        "C1Sf_nei",
        C1_nei * mesh.magSf()
    ),

     ap1
    (
      IOobject
      (
          "ap1",
          mesh.time().timeName(),
          //runTime.timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap1", dimensionSet(0, 3, -1, 0, 0, 0, 0), 1.0)
    ),

     am1
    (
      IOobject
      (
          "am1",
          mesh.time().timeName(),
          //runTime.timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap1", dimensionSet(0, 3, -1, 0, 0, 0, 0), -1.0)
    ),

     alpha1_own
    (
        "alpha1_own",
        ap1/(ap1 - am1)
    ),

     aSf1
    (
        "aSf1",
        am1*alpha1_own
    ),

     alpha1_nei
    (
        "alpha1_nei",
        1.0 - alpha1_own
    ),

     phi1_own
    (
        "phi1_own",
        phi_*rho1_own * 0.0
    ),

     phi1_nei
    (
        "phi1_nei",
        phi_*rho1_own * 0.0
    ),

    phi1_
    (
        "phi1",
        (phi1_own+phi1_nei)/(alpha1_own*rho1_own+alpha1_nei*rho1_nei)
    ),

/*******************************Region Two**********************************/

    rho2_own
    (
        "rho2own",
        fvc::interpolate(rho2_, own_, "reconstruct(rho2)")
    ),

    rho2_nei
    (
        "rho2nei",
        fvc::interpolate(rho2_, nei_, "reconstruct(rho2)")
    ),

    phi2v_own
    (
        "phi2v_own",
        (fvc::interpolate(rho2_*U_, own_, "reconstruct(U)") / rho2_own) & mesh.Sf()
    ),

    phi2v_nei
    (
        "phi2v_nei",
        (fvc::interpolate(rho2_*U_, nei_, "reconstruct(U)") / rho2_nei) & mesh.Sf()
    ),

     C2_own
    (
        "C2_own",
        fvc::interpolate(sqrt(C2_), own_, "reconstruct(phi)")
    ),

     C2_nei
    (
        "C2_nei",
        fvc::interpolate(sqrt(C2_), nei_, "reconstruct(phi)")
    ),

     C2Sf_own
    (
        "C2Sf_own",
        C2_own * mesh.magSf()
    ),

     C2Sf_nei
    (
        "C2Sf_nei",
        C2_nei * mesh.magSf()
    ),

     ap2
    (
      IOobject
      (
          "ap2",
          mesh.time().timeName(),
          //runTime.timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap2", dimensionSet(0, 3, -1, 0, 0, 0, 0), 1.0)
    ),

     am2
    (
      IOobject
      (
          "am2",
          mesh.time().timeName(),
          //runTime.timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap2", dimensionSet(0, 3, -1, 0, 0, 0, 0), -1.0)
    ),

     alpha2_own
    (
        "alpha2_own",
        ap2/(ap2 - am2)
    ),

     aSf2
    (
        "aSf2",
        am2*alpha2_own
    ),

     alpha2_nei
    (
        "alpha_nei",
        1.0 - alpha2_own
    ),

     phi2_own
    (
        "phi2_own",
        phi_*rho2_own * 0.0
    ),

     phi2_nei
    (
        "phi2_nei",
        phi_*rho2_own * 0.0
    ),

    phi2_
    (
        "phi2",
        (phi2_own+phi2_nei)/(alpha2_own*rho2_own+alpha2_nei*rho2_nei)
    ),

/***********************Tadmor-Kurganov Scheme*******************************/

/*************************Pressure Equation*********************************/

    psi1_own
    (
        "psi1_own",
        fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    psi1_nei
    (
        "psi1_nei",
        fvc::interpolate(psi1_, nei_, "reconstruct(psi)")
    ),

    psiU1_own
    (
        "psiU1_own",
//        fvc::interpolate(psi1_*HbyA_, own_, "reconstruct(U)")
         fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)")

    ),

    psiU1_nei
    (
        "psiU1_nei",
//        fvc::interpolate(psi1_*HbyA_, own_, "reconstruct(U)")
        fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)")
    ),

    aphi1v_own
    (
      "aphi1v_own",
      phi1v_own
    ),

    aphi1v_nei
    (
      "aphi1v_nei",
      phi1v_nei
    ),

    phi1d_own
    (
        "phi1d_own",
        aphi1v_own * psi1_own
    ),

    phi1d_nei
    (
        "phi1d_nei",
        aphi1v_nei * psi1_nei
    ),

    Dp1_own
    (
        "Dp1_own",
        alpha1_own*fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp1_nei
    (
        "Dp1_nei",
        alpha1_own*fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    psi2_own
    (
        "psi2_own",
        fvc::interpolate(psi2_, own_, "reconstruct(psi)")
    ),

    psi2_nei
    (
        "psi2_nei",
        fvc::interpolate(psi2_, nei_, "reconstruct(psi)")
    ),

    psiU2_own
    (
        "psiU2_own",
//        fvc::interpolate(psi2_*HbyA_, own_, "reconstruct(U)")
        fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)")
    ),

    psiU2_nei
    (
        "psiU2_nei",
//        fvc::interpolate(psi2_*HbyA_, own_, "reconstruct(U)")
         fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)")
    ),

    aphi2v_own
    (
      "aphi2v_own",
      phi2v_own
    ),

    aphi2v_nei
    (
      "aphi2v_nei",
      phi2v_nei
    ),

    phi2d_own
    (
        "phi2d_own",
        aphi2v_own * psi2_own
    ),

    phi2d_nei
    (
        "phi2d_nei",
        aphi2v_nei * psi2_nei
    ),

    Dp2_own
    (
        "Dp2_own",
        alpha2_own*fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp2_nei
    (
        "Dp2_nei",
        alpha2_own*fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    pEqn1_own
    (
        fvm::div(phi1d_own,p_) - fvm::laplacian(Dp1_own, p_)
    ),

    pEqn1_nei
    (
        fvm::div(phi1d_nei,p_) - fvm::laplacian(Dp1_nei, p_)
    ),

    pEqn2_own
    (
        fvm::div(phi2d_own,p_) - fvm::laplacian(Dp2_own, p_)
    ),

    pEqn2_nei
    (
        fvm::div(phi2d_nei,p_) - fvm::laplacian(Dp2_nei, p_)
    ),

    p_own
    (
        "p_own",
        fvc::interpolate(p_, own_, "reconstruct(p)")
    ),

    p_nei
    (
        "p_nei",
        fvc::interpolate(p_, own_, "reconstruct(p)")
    ),

    gradp1
    (
      "gradp1",
      fvc::div((alpha1_own*p_own + alpha1_nei*p_nei)*mesh.Sf())
    ),

    gradp2
    (
      "gradp2",
      fvc::div((alpha2_own*p_own + alpha2_nei*p_nei)*mesh.Sf())
    ),

/************************kappa Function*********************************/
/*
    kappaFuncPtr
    (
       fv::kappaFunction::New
       (
         "kappaFunction",
          mesh.solutionDict().subDict("PIMPLE").subDict("kappaFunction"),
          mesh
        )
    ),
*/
    kappa
    (
      IOobject
      (
          "kappa",
          mesh.time().constant(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      mesh,
      dimless
    ),

    amaxSf
    (
      "amaxSf",
       max(mag(aphi1v_own), mag(aphi1v_nei))
    ),

     uMagSf
     (
          IOobject
          (
            "uMagSf",
            mesh.time().constant(),
            mesh
          ),
        mesh.magSf() * 1.0
     ),

     divDevRhoReff1_
     (
       - fvm::laplacian(mu1_, U_)
       - fvc::div((mu1_)*dev2(Foam::T(fvc::grad(U_))))
     ),

     divDevRhoReff2_
     (
       - fvm::laplacian(mu2_, U_)
       - fvc::div((mu2_)*dev2(Foam::T(fvc::grad(U_))))
     ),

     Tviscosity1
     (
       - fvm::laplacian(alpha1_*Cp1_, T_)
       - 0.5*fvc::laplacian(alpha1_, magSqr(U_))
     ),

     Tviscosity2
     (
       - fvm::laplacian(alpha2_*Cp2_, T_)
       - 0.5*fvc::laplacian(alpha2_, magSqr(U_))
     ),

     devRhoReff1_
     (
       (-(alpha1_)*dev(twoSymm(fvc::grad(U_))))
     ),

     devRhoReff2_
     (
       (-(alpha2_)*dev(twoSymm(fvc::grad(U_))))
     ),

     TSource1_
     (
       fvc::div((linearInterpolate((-devRhoReff1_) & U_) & mesh.Sf())())
     ),

     TSource2_
     (
       fvc::div((linearInterpolate((-devRhoReff2_) & U_) & mesh.Sf())())
     ),

     psi1f
     (
       "psi1f",
       linearInterpolate(psi1_)
     ),

     psi2f
     (
       "psi2f",
       linearInterpolate(psi1_)
     ),

    rho01_own
    (
        fvc::interpolate(rho01_, own_, "reconstruct(rho2)")
    ),

    rho01_nei
    (
        fvc::interpolate(rho01_, nei_, "reconstruct(rho2)")
    ),

    rho02_own
    (
        fvc::interpolate(rho02_, own_, "reconstruct(rho2)")
    ),

    rho02_nei
    (
        fvc::interpolate(rho02_, nei_, "reconstruct(rho2)")
    ),

     phi01d_own
     (
       aphi1v_own * rho01_own
     ),

     phi01d_nei
     (
       aphi1v_nei * rho01_nei
     ),

     phi02d_own
     (
       aphi2v_own * rho02_own
     ),

     phi02d_nei
     (
       aphi2v_nei * rho02_nei
     ),

     scAlpha_
     (
       dimensioned< scalar >("scAlpha", *this)
     ),

     phic_
     (
       mag(phi_/mesh.magSf())
     ),

     E_
     (
       (
         fvc::ddt(rho_)
         +
         fvc::div(phi1_own+phi1_nei,volumeFraction1_)
         +
         fvc::div(phi1_own+phi2_nei,volumeFraction2_)
       )
     ),

     E1_
     (
       (
         fvc::ddt(rho1_)
         +
         fvc::div(phi1_own+phi1_nei)
       )
     ),

     E2_
     (
       (
         fvc::ddt(rho2_)
         +
         fvc::div(phi2_own+phi2_nei)
       )
     ),

     EF_
     (
        volumeFraction1_*E1_ + volumeFraction2_*E2_
     ),

     Ep_
     (
       0*EF_
     ),

     ddtvF1_
     (
       fvc::ddt(volumeFraction1_)
     ),

     ddtvF2_
     (
       fvc::ddt(volumeFraction2_)
     ),

     vF1face_
     (
       fvc::interpolate(volumeFraction1_, "reconstruct(volumeFraction1)")
     ),

     vF2face_
     (
       1.0 - vF1face_
     ),

     vFPhi1_
     (
       phi_*vF1face_
     ),

     vFPhi2_
     (
       phi_ - vFPhi1_
     ),

     divU_
     (
       volumeFraction1_*E1_/rho1_
       +
       volumeFraction2_*E2_/rho2_
       -
       ddtvF1_ - fvc::div(vFPhi1_)
       -
       ddtvF2_ - fvc::div(vFPhi2_)
     )


{
    Info<< "\nConstructor is working\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTwoPhaseCentralFoam::~interTwoPhaseCentralFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interTwoPhaseCentralFoam::saveOld()
    {

      volumeFraction1_.oldTime();
      rho1_.oldTime();
      rho2_.oldTime();
      U_.oldTime();
      T_.oldTime();
      p_.oldTime();
      psi1_.oldTime();
      psi2_.oldTime();
      rho01_.oldTime();
      rho02_.oldTime();
      rho1AmpU_.oldTime();
      rho2AmpU_.oldTime();
      TSource1_.oldTime();
      TSource2_.oldTime();
      ddtvF1_.oldTime();
      ddtvF2_.oldTime();

    }



// * * * * * * * * * * * * * * * Solving Functions  * * * * * * * * * * * * * //

    void Foam::interTwoPhaseCentralFoam::alpha1Eqnsolve()
    {

      vF1face_ = fvc::interpolate(volumeFraction1_, "reconstruct(volumeFraction1)");
      vF2face_ = 1.0 - vF1face_;
      vFPhi1_ = phi_*vF1face_;
      vFPhi2_ = phi_ - vFPhi1_;

      fvScalarMatrix alpha1Eqn
      (
/*
        fvm::ddt(volumeFraction1_)
        +
        fvc::div(vFPhi1_)
        ==
        (1+K_)*fvc::div(phi_)*volumeFraction1_
*/
        fvm::ddt(volumeFraction1_)
        +
        fvc::div(phi_, volumeFraction1_)
        ==
        (1+K_)*fvc::div(phi_)*volumeFraction1_

      );

      alpha1Eqn.solve();

      volumeFraction2_ = 1 - volumeFraction1_;



      ddtvF1_ = fvc::ddt(volumeFraction1_);
      ddtvF2_ = fvc::ddt(volumeFraction2_);

      Info<< "max: volumeFraction1 " << max(volumeFraction1_).value()
          << " min: " << min(volumeFraction1_).value()
          << nl << endl;
    }


    void Foam::interTwoPhaseCentralFoam::UEqnsolve()
    {

      pressureGradient();
      divDevRhoReff();
      fvVectorMatrix UEqn
      (
       volumeFraction1_ *
        (fvm::ddt(rho1_,U_) + fvm::div(phi1_own,U_)
        +
        fvm::div(phi1_nei,U_) + gradp1 + divDevRhoReff1_)
       +
       volumeFraction2_ *
       (fvm::ddt(rho2_,U_) + fvm::div(phi2_own,U_)
       +
       fvm::div(phi2_nei,U_) + gradp2 + divDevRhoReff2_)
       );

    UEqn.solve();

    rbyA_  = 1.0 / UEqn.A();
    HbyA_ = UEqn.H() * rbyA_;
    HbyA_.boundaryFieldRef() == U_.boundaryField();

    }


    void Foam::interTwoPhaseCentralFoam::TEqnsolve()
    {

      fvScalarMatrix TEqn
      (
        volumeFraction1_ *
        (
          Cp1_*fvm::ddt(rho1_,T_)
//          +
//          0.5*fvc::ddt(rho1AmpU_)
          +
          Cp1_*fvm::div(phi1_own,T_) + Cp1_*fvm::div(phi1_nei,T_)
          +
          0.5*(fvc::div(phi1_own,magSqr(U_)) + fvc::div(phi1_nei,magSqr(U_)))
          +
          Tviscosity1
          -
          TSource1_
        )
        +
        volumeFraction2_ *
        (
           Cp2_*fvm::ddt(rho2_,T_)
//           +
//           0.5*fvc::ddt(rho2AmpU_)
          +
          Cp2_*fvm::div(phi2_own,T_) + Cp2_*fvm::div(phi2_nei,T_)
          +
          0.5*(fvc::div(phi2_own,magSqr(U_)) + fvc::div(phi2_nei,magSqr(U_)))
          +
          Tviscosity2
          -
          TSource2_
        )
//        ==
//        fvc::ddt(p_)

      );

      TEqn.solve();

    }


    void Foam::interTwoPhaseCentralFoam::rho1Eqnsolve()
    {
      fvScalarMatrix rho1Eqn
      (
        (
        fvm::ddt(rho1_)
        +
        fvc::div(phi1_own) + fvc::div(phi1_nei)
        )
      );

      rho1Eqn.solve();

    }

    void Foam::interTwoPhaseCentralFoam::rho2Eqnsolve()
    {
      fvScalarMatrix rho2Eqn
      (
        (
        fvm::ddt(rho2_)
        +
        fvc::div(phi2_own) + fvc::div(phi2_nei)
        )
      );

      rho2Eqn.solve();

    }

    void Foam::interTwoPhaseCentralFoam::pEqnsolve()
    {

      pEqn1_own = fvc::div(phi01d_own) + fvm::div(phi1d_own,p_) - fvm::laplacian(Dp1_own, p_);
      pEqn1_nei = fvc::div(phi01d_nei) + fvm::div(phi1d_nei,p_) - fvm::laplacian(Dp1_nei, p_);

      pEqn2_own = fvc::div(phi02d_own) + fvm::div(phi2d_own,p_) - fvm::laplacian(Dp2_own, p_);
      pEqn2_nei = fvc::div(phi02d_nei) + fvm::div(phi2d_nei,p_) - fvm::laplacian(Dp2_nei, p_);

        fvScalarMatrix pEqn
        (

            volumeFraction1_*
            (fvm::ddt(psi1_,p_) + fvc::ddt(rho01_)
            + pEqn1_own + pEqn1_nei
            )
            +
            volumeFraction2_*
            (fvm::ddt(psi2_,p_) + fvc::ddt(rho02_)
            + pEqn2_own + pEqn2_nei
            )
        );

        pEqn.solve();
/*
        Ep_ =volumeFraction2_ *( (pEqn2_own & p_)
        + (pEqn2_nei & p_) + fvc::ddt(psi2_,p_) + fvc::ddt(rho02_));
        Ep_.rename("Ep");
        Ep_.write();
*/
    }

// * * * * * * * * * * * * * * * Mid-Step Functions * * * * * * * * * * * * * //

    void Foam::interTwoPhaseCentralFoam::Update()
    {

      C1_ = gamma1_*R1_*T_;
      C2_ = gamma2_*R2_*T_;

      Z1_ = rho1_*C1_;
      Z2_ = rho2_*C2_;
      K_ = (volumeFraction1_*volumeFraction2_*(Z2_-Z1_))
      /
      (Z1_*volumeFraction2_+Z2_*volumeFraction1_);

    }

    void Foam::interTwoPhaseCentralFoam::UpdateVelocity()
    {

      rho1AmpU_ = magSqr(U_) * rho1_;
      rho2AmpU_ = magSqr(U_) * rho2_;

    }


    void Foam::interTwoPhaseCentralFoam::rhoU1TKS_()
    {

      rho1_own     = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei     = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      phi1v_own    = (fvc::interpolate(rho1_*U_, own_, "reconstruct(U)") / rho1_own) & U_.mesh().Sf();
      phi1v_nei    = (fvc::interpolate(rho1_*U_, nei_, "reconstruct(U)") / rho1_nei) & U_.mesh().Sf();

      C1_own      = fvc::interpolate(sqrt(C1_), own_, "reconstruct(psi1)");
      C1_nei      = fvc::interpolate(sqrt(C1_), nei_, "reconstruct(psi1)");

      C1Sf_own = C1_own * U_.mesh().magSf();
      C1Sf_own.setOriented(true);
      C1Sf_nei = C1_nei * U_.mesh().magSf();
      C1Sf_nei.setOriented(true);

      ap1 = max(max(phi1v_own + C1Sf_own, phi1v_nei + C1Sf_nei), v_zero);
      am1 = min(min(phi1v_own - C1Sf_own, phi1v_nei - C1Sf_nei), v_zero);


      alpha1_own   = ap1/(ap1 - am1);
      aSf1         = am1*alpha1_own;
      alpha1_nei   = 1.0 - alpha1_own;

      phi1_own = rho1_own * (alpha1_own*phi1v_own - aSf1);
      phi1_nei = rho1_nei * (alpha1_nei*phi1v_nei + aSf1);


    }

    void Foam::interTwoPhaseCentralFoam::rhoU2TKS_()
    {

      rho2_own     = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
      rho2_nei     = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

      phi2v_own    = (fvc::interpolate(rho2_*U_, own_, "reconstruct(U)") / rho2_own) & U_.mesh().Sf();
      phi2v_nei    = (fvc::interpolate(rho2_*U_, nei_, "reconstruct(U)") / rho2_nei) & U_.mesh().Sf();

      C2_own      = fvc::interpolate(sqrt(C2_), own_, "reconstruct(psi2)");
      C2_nei      = fvc::interpolate(sqrt(C2_), nei_, "reconstruct(psi2)");

      C2Sf_own = C2_own * U_.mesh().magSf();
      C2Sf_own.setOriented(true);
      C2Sf_nei = C2_nei * U_.mesh().magSf();
      C2Sf_nei.setOriented(true);

      ap2 = max(max(phi2v_own + C2Sf_own, phi2v_nei + C2Sf_nei), v_zero);
      am2 = min(min(phi2v_own - C2Sf_own, phi2v_nei - C2Sf_nei), v_zero);


      alpha2_own   = ap2/(ap2 - am2);
      aSf2         = am2*alpha2_own;
      alpha2_nei   = 1.0 - alpha2_own;

      phi2_own = rho2_own * (alpha2_own*phi2v_own - aSf2);
      phi2_nei = rho2_nei * (alpha2_nei*phi2v_nei + aSf2);

    }


    void Foam::interTwoPhaseCentralFoam::Flux()
    {

      phi1_own = pEqn1_own.flux() + phi01d_own;
      phi1_nei = pEqn1_nei.flux() + phi01d_nei;

      phi2_own = pEqn2_own.flux() + phi02d_own;
      phi2_nei = pEqn2_nei.flux() + phi02d_nei;

    }

    void Foam::interTwoPhaseCentralFoam::alphaPhi()
    {

        phi_ =
        (phi1_own+phi1_nei)/(alpha1_own*rho1_own+alpha1_nei*rho1_nei)
//        +
//        (phi2_own+phi2_nei)/(alpha2_own*rho2_own+alpha2_nei*rho2_nei)
        ;

/*
        surfaceScalarField rbyAf = fvc::interpolate(rbyA_);

        divU();
        Foam::CorrectPhi
        (
          U_,
          phi_,
          p_,
          rbyAf,
          divU_,
          pimple_
        );
*/

        Info << "C1: " << max(sqrt(C1_)).value() << endl;
        Info << "C2: " << max(sqrt(C2_)).value() << endl;
    }


    void Foam::interTwoPhaseCentralFoam::Velocity()
    {

      pressureGradient();
//      U_ = HbyA_-volumeFraction1_*rbyA_*gradp1-volumeFraction2_*rbyA_*gradp2;
      U_ = HbyA_-rbyA_*gradp2;
      U_.correctBoundaryConditions();

    }

    void Foam::interTwoPhaseCentralFoam::DensityThermo()
    {
        rho1_ = psi1_*p_+rho01Sc_;
        rho2_ = psi2_*p_+rho02Sc_;
        rho1_.correctBoundaryConditions();
        rho2_.correctBoundaryConditions();

        p_own = fvc::interpolate(p_, own_, "reconstruct(p)");
        p_nei = fvc::interpolate(p_, nei_, "reconstruct(p)");

    }

    void Foam::interTwoPhaseCentralFoam::DensityNull()
    {
      rho01_ = rho1_ - psi1_*p_;
      rho02_ = rho2_ - psi2_*p_;
    }

    void Foam::interTwoPhaseCentralFoam::DensityCorrection()
    {
      rho1_ = psi1_*p_ + rho01_;
      rho2_ = psi2_*p_ + rho02_;
    }

    void Foam::interTwoPhaseCentralFoam::Density()
    {
      rho_ = volumeFraction1_ * rho1_ + volumeFraction2_ * rho2_;
    }

    void Foam::interTwoPhaseCentralFoam::ThermoCoefficient()
    {
      psi1_ = 1 / (R1_ * T_);
      psi2_ = 1 / (R2_ * T_);

    }

    void Foam::interTwoPhaseCentralFoam::Initialize()
    {
      psi1_ = (1 / (R1_ * T_));
      psi2_ = (1 / (R2_ * T_));

      rho01_ = (rho01Sc_*rho1_)/rho1_;
      rho02_ = (rho02Sc_*rho2_)/rho2_;

      rho1_ = psi1_*p_ + rho01_;
      rho2_ = psi2_*p_ + rho02_;
      rho1_.correctBoundaryConditions();
      rho2_.correctBoundaryConditions();
      rho_ = volumeFraction1_ * rho1_ + volumeFraction2_ * rho2_;

      gamma1_ = Cp1_ / (Cp1_ - (R_/molM1_));
      gamma2_ = Cp2_ / (Cp2_ - (R_/molM2_));

      Update();

      alpha();
      UpdateCentralWeights();
      UpdateKappa();
      divDevRhoReff();

      UEqn();

      UpdateCentralFields();

      pressureGradient();

      devRhoReff();

      rho1AmpU_ = magSqr(U_) * rho1_;
      rho2AmpU_ = magSqr(U_) * rho2_;

      Tviscosity1 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_)
                - 0*fvc::laplacian(alpha1_, magSqr(U_)));
      Tviscosity2 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_)
                - 0*fvc::laplacian(alpha1_, magSqr(U_)));

      TSource1_
      = 0*
      fvc::div
      ((linearInterpolate((-devRhoReff1_) & U_) & U_.mesh().Sf())());
      TSource2_ = 0*
      fvc::div
      ((linearInterpolate((-devRhoReff2_) & U_) & U_.mesh().Sf())());

    }

    void Foam::interTwoPhaseCentralFoam::pressureGradient()
    {

      p_own = fvc::interpolate(p_, own_, "reconstruct(p)");
      p_nei = fvc::interpolate(p_, nei_, "reconstruct(p)");

      gradp1 = fvc::div((alpha1_own*p_own + alpha1_nei*p_nei)*U_.mesh().Sf());
      gradp2 = fvc::div((alpha2_own*p_own + alpha2_nei*p_nei)*U_.mesh().Sf());

    }

    void Foam::interTwoPhaseCentralFoam::CallFlux()
    {
      rhoU1TKS_();
      rhoU2TKS_();
    }

    void Foam::interTwoPhaseCentralFoam::UEqn()
    {
      divDevRhoReff();    //Add a viscous term to UEqn
      pressureGradient();

      fvVectorMatrix UEqn
      (

       volumeFraction1_ *
       (fvm::ddt(rho1_,U_) + fvm::div(phi1_own,U_)
        +
       fvm::div(phi1_nei,U_) + divDevRhoReff1_)
       +
       volumeFraction2_ *
       (fvm::ddt(rho2_,U_) + fvm::div(phi2_own,U_)
        + fvm::div(phi2_nei,U_) + divDevRhoReff2_)

       );

      rbyA_  = 1.0 / UEqn.A();
      HbyA_ = UEqn.H() * rbyA_;
      HbyA_.boundaryFieldRef() == U_.boundaryField();
    }

    void Foam::interTwoPhaseCentralFoam::UpdateCentralWeights()
    {

      rho1_own     = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei     = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      rho2_own     = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
      rho2_nei     = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

      phi1v_own    = (fvc::interpolate(U_, own_, "reconstruct(U)")) & U_.mesh().Sf();
      phi1v_nei    = (fvc::interpolate(U_, nei_, "reconstruct(U)")) & U_.mesh().Sf();

      phi2v_own    = (fvc::interpolate(U_, own_, "reconstruct(U)")) & U_.mesh().Sf();
      phi2v_nei    = (fvc::interpolate(U_, nei_, "reconstruct(U)")) & U_.mesh().Sf();
      /*

      rho1_own     = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei     = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      rho2_own     = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
      rho2_nei     = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

      phi1v_own    = (fvc::interpolate(rho1_*U_, own_, "reconstruct(U)") / rho1_own) & U_.mesh().Sf();
      phi1v_nei    = (fvc::interpolate(rho1_*U_, nei_, "reconstruct(U)") / rho1_nei) & U_.mesh().Sf();

      phi2v_own    = (fvc::interpolate(rho2_*U_, own_, "reconstruct(U)") / rho2_own) & U_.mesh().Sf();
      phi2v_nei    = (fvc::interpolate(rho2_*U_, nei_, "reconstruct(U)") / rho2_nei) & U_.mesh().Sf();


      C2_own      = fvc::interpolate(sqrt(C2_), own_, "reconstruct(psi2)");
      C2_nei      = fvc::interpolate(sqrt(C2_), nei_, "reconstruct(psi2)");

      C1_own      = fvc::interpolate(sqrt(C1_), own_, "reconstruct(psi1)");
      C1_nei      = fvc::interpolate(sqrt(C1_), nei_, "reconstruct(psi1)");
*/
      speedOfSound();
      C2_own      = fvc::interpolate((C_), own_, "reconstruct(psi2)");
      C2_nei      = fvc::interpolate((C_), nei_, "reconstruct(psi2)");

      C1_own      = fvc::interpolate((C_), own_, "reconstruct(psi1)");
      C1_nei      = fvc::interpolate((C_), nei_, "reconstruct(psi1)");


      C1Sf_own = C1_own * U_.mesh().magSf();
      C1Sf_own.setOriented(true);
      C1Sf_nei = C1_nei * U_.mesh().magSf();
      C1Sf_nei.setOriented(true);


      ap1 = max(max(phi1v_own + C1Sf_own, phi1v_nei + C1Sf_nei), v_zero);
      am1 = min(min(phi1v_own - C1Sf_own, phi1v_nei - C1Sf_nei), v_zero);


      C2Sf_own = C2_own * U_.mesh().magSf();
      C2Sf_own.setOriented(true);
      C2Sf_nei = C2_nei * U_.mesh().magSf();
      C2Sf_nei.setOriented(true);

      ap2 = max(max(phi2v_own + C2Sf_own, phi2v_nei + C2Sf_nei), v_zero);
      am2 = min(min(phi2v_own - C2Sf_own, phi2v_nei - C2Sf_nei), v_zero);

      alpha1_own   = ap1/(ap1 - am1);
      aSf1         = am1*alpha1_own;
      alpha1_nei   = 1.0 - alpha1_own;

      alpha2_own   = ap2/(ap2 - am2);
      aSf2         = am2*alpha2_own;
      alpha2_nei   = 1.0 - alpha2_own;

    }

    void Foam::interTwoPhaseCentralFoam::UpdateCentralFields()
    {

      rho1_own     = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei     = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      rho2_own     = fvc::interpolate(rho2_, own_, "reconstruct(rho1)");
      rho2_nei     = fvc::interpolate(rho2_, nei_, "reconstruct(rho1)");


      psi1_own = fvc::interpolate(psi1_, own_, "reconstruct(psi)");
      psi1_nei = fvc::interpolate(psi1_, nei_, "reconstruct(psi)");

      psi2_own = fvc::interpolate(psi2_, own_, "reconstruct(psi)");
      psi2_nei = fvc::interpolate(psi2_, nei_, "reconstruct(psi)");

/*
      psi1f = linearInterpolate(psi1_);
      psi2f = linearInterpolate(psi2_);

      psi1_own     = fvc::interpolate(psi1_, own_, "reconstruct(psi)")*kappa
                    + (1.0 - kappa)*psi1f;
      psi1_nei     = fvc::interpolate(psi1_, nei_, "reconstruct(psi)")*kappa
                    + (1.0 - kappa)*psi1f;

      psi2_own     = fvc::interpolate(psi2_, own_, "reconstruct(psi)")*kappa
                    + (1.0 - kappa)*psi2f;
      psi2_nei     = fvc::interpolate(psi2_, nei_, "reconstruct(psi)")*kappa
                    + (1.0 - kappa)*psi2f;
*/
/*
      psiU1_own    = fvc::interpolate(psi1_*HbyA_, own_, "reconstruct(U)");
      psiU1_nei    = fvc::interpolate(psi1_*HbyA_, nei_, "reconstruct(U)");

      psiU2_own    = fvc::interpolate(psi2_*HbyA_, own_, "reconstruct(U)");
      psiU2_nei    = fvc::interpolate(psi2_*HbyA_, nei_, "reconstruct(U)");

      phi1v_own    = (psiU1_own / psi1_own) & U_.mesh().Sf();
      phi1v_nei    = (psiU1_nei / psi1_nei) & U_.mesh().Sf();

      phi2v_own    = (psiU2_own / psi2_own) & U_.mesh().Sf();
      phi2v_nei    = (psiU2_nei / psi2_nei) & U_.mesh().Sf();
*/
      psiU1_own    = fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)");
      psiU1_nei    = fvc::interpolate(rho1_*HbyA_, nei_, "reconstruct(U)");

      psiU2_own    = fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)");
      psiU2_nei    = fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)");

      phi1v_own    = (psiU1_own / rho1_own) & U_.mesh().Sf();
      phi1v_nei    = (psiU1_nei / rho1_nei) & U_.mesh().Sf();

      phi2v_own    = (psiU2_own / rho2_own) & U_.mesh().Sf();
      phi2v_nei    = (psiU2_nei / rho2_nei) & U_.mesh().Sf();

      aphi1v_own   = alpha1_own*phi1v_own - aSf1;
      aphi1v_nei   = alpha1_nei*phi1v_nei + aSf1;

      aphi2v_own   = alpha2_own*phi2v_own - aSf2;
      aphi2v_nei   = alpha2_nei*phi2v_nei + aSf2;

      phi1d_own    = aphi1v_own * psi1_own;
      phi1d_nei    = aphi1v_nei * psi1_nei;

      phi2d_own    = aphi2v_own * psi2_own;
      phi2d_nei    = aphi2v_nei * psi2_nei;

      Dp1_own = alpha1_own*fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)");
      Dp1_nei = alpha1_nei*fvc::interpolate(rho1_*rbyA_, nei_, "reconstruct(Dp)");

      Dp2_own = alpha2_own*fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)");
      Dp2_nei = alpha2_nei*fvc::interpolate(rho2_*rbyA_, nei_, "reconstruct(Dp)");
/*
      rho01_own     = fvc::interpolate(rho01_, own_, "reconstruct(rho1)");
      rho01_nei     = fvc::interpolate(rho01_, nei_, "reconstruct(rho1)");

      rho02_own     = fvc::interpolate(rho02_, own_, "reconstruct(rho1)");
      rho02_nei     = fvc::interpolate(rho02_, nei_, "reconstruct(rho1)");
*/
      phi01d_own    = aphi1v_own * rho01Sc_;
      phi01d_nei    = aphi1v_nei * rho01Sc_;

      phi02d_own    = aphi2v_own * rho02Sc_;
      phi02d_nei    = aphi2v_nei * rho02Sc_;

/*
      phi01d_own    = alpha1_own*phi1v_own * rho01Sc_;
      phi01d_nei    = alpha1_nei*phi1v_nei * rho01Sc_;

      phi02d_own    = alpha2_own*phi2v_own * rho02Sc_;
      phi02d_nei    = alpha2_nei*phi2v_nei * rho02Sc_;
*/

    }


    void Foam::interTwoPhaseCentralFoam::UpdateKappa()
    {/*
      uMagSf = U_.mesh().magSf();
      phi1_ = phi1_own + phi1_nei;

      aphi1v_own = phi1_own / rho1_own;
      aphi1v_nei = phi1_nei / rho1_nei;
      dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(U_.mesh().magSf()).value());
//      amaxSmall = SMALL * min(U_.mesh().magSf()).value();
      amaxSf = max(mag(aphi1v_own), mag(aphi1v_nei)) + amaxSmall;
      amaxSf.setOriented(true);

      kappaFuncPtr().update();
      kappa = kappaFuncPtr().kappa()();

      Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;

//      phi1_own += (1.0 - kappa) * phi1_nei;
//      phi1_nei = kappa * phi1_nei;

//      phi2_own += (1.0 - kappa) * phi2_nei;
//      phi2_nei = kappa * phi2_nei;

        phi1_own += 0*phi1_nei;
        phi1_nei = 1*phi1_nei;

        phi2_own += 0*phi2_nei;
        phi2_nei = 1*phi2_nei;
*/
    }

    void Foam::interTwoPhaseCentralFoam::divDevRhoReff()
    {

        divDevRhoReff1_ =
        (
        - fvm::laplacian(mu1_, U_)
        - fvc::div((mu1_)*dev2(Foam::T(fvc::grad(U_))))
        );

        divDevRhoReff2_ =
        (
        - fvm::laplacian(mu2_, U_)
        - fvc::div((mu2_)*dev2(Foam::T(fvc::grad(U_))))
        );

    }


    void Foam::interTwoPhaseCentralFoam::Tviscosity()
    {
        alpha();
        Tviscosity1 =
        (
          - fvm::laplacian(alpha1_*Cp1_, T_)
          - 0.5*fvc::laplacian(alpha1_, magSqr(U_))
        );

        Tviscosity2 =
        (
          - fvm::laplacian(alpha2_*Cp2_, T_)
          - 0.5*fvc::laplacian(alpha2_, magSqr(U_))
        );

    }

    void Foam::interTwoPhaseCentralFoam::alpha()
    {
      alpha1_ = mu1_/Pr_;
      alpha2_ = mu2_/Pr_;
      alpha_ = volumeFraction1_*alpha1_+volumeFraction2_*alpha2_;
    }

    void Foam::interTwoPhaseCentralFoam::devRhoReff()
    {
      alpha();
      devRhoReff1_ = (-(alpha1_)*dev(twoSymm(fvc::grad(U_))));
      devRhoReff2_ = (-(alpha2_)*dev(twoSymm(fvc::grad(U_))));

    }

    void Foam::interTwoPhaseCentralFoam::TViscousitySource()
    {

      rho1AmpU_ = magSqr(U_) * rho1_;
      rho2AmpU_ = magSqr(U_) * rho2_;

      devRhoReff();

      TSource1_
      =
      fvc::div
      ((linearInterpolate((-devRhoReff1_) & U_) & U_.mesh().Sf())())
      +
       fvc::ddt(p_)
      - 0.5*fvc::ddt(rho1AmpU_);

      TSource2_
      =
      fvc::div
      ((linearInterpolate((-devRhoReff2_) & U_) & U_.mesh().Sf())())
      +
       fvc::ddt(p_)
      - 0.5*fvc::ddt(rho2AmpU_);
    }

    void Foam::interTwoPhaseCentralFoam::info()
    {

      E_ =
      (
        fvc::ddt(rho_)
        +
        fvc::div(phi1_own+phi1_nei,volumeFraction1_)
        +
        fvc::div(phi1_own+phi2_nei,volumeFraction2_)
      );

      E1_ =
      (
        fvc::ddt(rho01_)+fvc::ddt(psi1_,p_)
        +
        fvc::div(phi1_own+phi1_nei)
      );

      E2_ =
      (
        fvc::ddt(rho2_)
        +
        fvc::div(phi2_own+phi2_nei)
      );

      EF_ = volumeFraction1_*E1_ + volumeFraction2_*E2_;

      EF_.rename("EF");
      EF_.write();

      Info<< " max: E " << max(E_).value()
          << " min: " << min(E_).value()
          << " max: E1 " << max(E1_).value()
          << " min: " << min(E1_).value()
          << " max: E2 " << max(E2_).value()
          << " min: " << min(E2_).value()
          << " max: EF " << max(EF_).value()
          << " min: " << min(EF_).value()
          << nl << endl;


    }

    void Foam::interTwoPhaseCentralFoam::alphaPhiCorrect()
    {
      phic_ = mag(phi_/U_.mesh().magSf());
    }

    void Foam::interTwoPhaseCentralFoam::divU()
    {
      E1_ =
      (
        fvc::ddt(rho01_)+fvc::ddt(psi1_,p_)
        +
        fvc::div(phi1_own+phi1_nei)
      );

      E2_ =
      (
        fvc::ddt(rho02_)+fvc::ddt(psi2_,p_)
        +
        fvc::div(phi2_own+phi2_nei)
      );

      divU_ =
        volumeFraction1_*E1_/rho1_
        +
        volumeFraction2_*E2_/rho2_
        -
        ddtvF1_ - fvc::div(vFPhi1_)
        -
        ddtvF2_ - fvc::div(vFPhi2_);

    }

    void Foam::interTwoPhaseCentralFoam::speedOfSound()
    {

      volScalarField rbypsiM = 1/(volumeFraction1_*psi1_+volumeFraction2_*psi2_);
      volScalarField psiM = 1/rbypsiM;
      volScalarField y1 = volumeFraction1_*(rho1_/rho_);
      volScalarField y2 = volumeFraction2_*(rho2_/rho_);
      volScalarField CpM = y1*Cp1_ + y2*Cp2_;
      volScalarField CvM = y1*(Cp1_/gamma1_) + y2*(Cp2_/gamma2_);
      volScalarField gammaM = CpM/CvM;
      C_ = sqrt(gammaM/psiM);

    }
