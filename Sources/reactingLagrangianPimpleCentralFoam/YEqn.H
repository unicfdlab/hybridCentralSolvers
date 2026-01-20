// mesh addressing might change, therefore, diffusive fluxes must
// be allocated again

forAll(Y, iCmpt)
{
    diffusiveFlux.set(iCmpt, phi*0.0);
}
#include "../reactingPimpleCentralFoam/YEqn.H"

