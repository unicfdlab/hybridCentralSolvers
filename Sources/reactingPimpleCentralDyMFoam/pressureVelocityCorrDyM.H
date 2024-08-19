if (updateEnergyInPISO) //update each iteration before pressure
{
    #include "YEqn.H"
    #include "hEqn.H"
}
#include "pEqn.H"
if (updateEnergyInPISO)
{
    #define PISOCENTRALFOAM_LTS
    if (mesh.moving())
    {
        mphi_own = alpha_own * rho_own * mesh_phi;
        mphi_nei = alpha_nei * rho_nei * mesh_phi;
        
        phi_own += mphi_own;
        phi_nei += mphi_nei;
        phi = phi_own + phi_nei;
    }

    // --- update blending function
    #include "updateKappa.H"
    // --- update mechanical fields
    #include "updateMechanicalFields.H"
}

