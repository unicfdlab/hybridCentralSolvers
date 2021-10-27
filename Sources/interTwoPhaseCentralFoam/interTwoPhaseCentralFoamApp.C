#include "fvCFD.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "interTwoPhaseCentralFoam.H"

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  scalar CoNum = 0.0;
  scalar meanCoNum = 0.0;

  pimpleControl pimple(mesh);
//  Info << "nCorrPIMPLE: " << pimple.nCorrPIMPLE()  << endl;

  interTwoPhaseCentralFoam Veronika(mesh);


  #include "createTimeControls.H"
  #include "setInitialDeltaT.H"

  Info<< "\nStarting time loop\n" << endl;

  Veronika.Initialize();

  while (runTime.run())
  {

/*
    scalarField sumC_
    (
        mag(Veronika.C())
    );
*/
    scalarField sumPhi_
    (
        fvc::surfaceSum(mag(Veronika.phiC()))().primitiveField()
    );

    CoNum = gMax((sumPhi_)/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        (gSum(sumPhi_)/gSum(mesh.V().field()))*runTime.deltaTValue();
        Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    #include "readTimeControls.H"

    #include "setDeltaT.H"

    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;
    Veronika.saveOld();
//    Veronika.alphaPhi();
//    Veronika.alphaEqnsolve();

//    Veronika.rho1Eqnsolve();
//    Veronika.rho2Eqnsolve();

    while (pimple.loop())
    {


      Veronika.alpha1Eqnsolve();
//      Veronika.alpha2Eqnsolve();

      Veronika.rho1Eqnsolve();
      Veronika.rho2Eqnsolve();

      Veronika.UEqn();          //Generate fvmatrix for UEqn (note: without grad(p_))

      Veronika.TEqnsolve();     //Solve TEqn

      Veronika.ThermoCoefficient(); //Update psi1_ = molM1_/(R_ * T_) and psi2_
      Veronika.Update();         //Calculate values of C1_, C2_, Z1_, Z2_, K_, and phi_
      Veronika.DensityThermo();

      Veronika.UpdateCentralWeights(); //Calculate fluxes (phi1_own and phi1_nei)
      Veronika.UpdateCentralFields();  //Calculate coefficients of pEqn: phi1d_own, phi1_nei, Dp1_own, and Dp2_nei

      Veronika.pEqnsolve();            //Solve pEqn
      Veronika.Flux();

      Veronika.DensityThermo();        //Update rho1_ and rho2_ through rho_ = psi_*p_
      Veronika.Density();
      Veronika.Velocity();

      Veronika.alphaPhi();

      Veronika.Tviscosity();
      Veronika.TViscousitySource();
    }

     runTime.write();



     Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;


  }

    Info<< "End\n" << endl;

    return 0;
}
