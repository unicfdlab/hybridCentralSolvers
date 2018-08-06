    if (updateEnergyInPISO) //update each iteration before pressure
    {
        #include "hEqn.H"
    }
    #include "pEqn.H"
    if (updateEnergyInPISO)
    {
        #define PISOCENTRALFOAM_LTS
        // --- update blending function
        #include "updateKappa.H"
        // --- update mechanical fields
        #include "updateMechanicalFields.H"
    }

