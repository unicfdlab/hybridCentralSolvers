#ifndef Two_Phase_Mixture_Update_Sonic_Speed_On_Faces_H
#define Two_Phase_Mixture_Update_Sonic_Speed_On_Faces_H

#include "updateSonicSpeedOnFaces.H"

template<>
void UpdateSonicSpeedOnFaces
(
    compressibleTwoPhaseMixtureThermo &thermo,
    const surfaceScalarField &own,
    const surfaceScalarField &nei,
    const surfaceScalarField &p_own,
    const surfaceScalarField &p_nei,
    const surfaceScalarField &phi,
    surfaceScalarField &cf_own,
    surfaceScalarField &cf_nei
)
{
    surfaceScalarField rho_liq_own = fvc::interpolate
    (
        thermo.thermoLiq().rho(), own, "reconstruct(rho)"
    );
    surfaceScalarField rho_liq_nei = fvc::interpolate
    (
        thermo.thermoLiq().rho(), nei, "reconstruct(rho)"
    );
    surfaceScalarField rho_gas_own = fvc::interpolate
    (
        thermo.thermoGas().rho(), own, "reconstruct(rho)"
    );
    surfaceScalarField rho_gas_nei = fvc::interpolate
    (
        thermo.thermoGas().rho(), nei, "reconstruct(rho)"
    );
    surfaceScalarField YbarLiqf    = fvc::interpolate
    (
        thermo.YbarLiq(), phi, "reconstruct(YbarLiq)"
    );

    const volScalarField& T = thermo.T();
    surfaceScalarField T_own =
        fvc::interpolate(T, own, "reconstruct(psi)");
    surfaceScalarField T_nei =
        fvc::interpolate(T, nei, "reconstruct(psi)");

    cf_own = thermo.Cfm(p_own, T_own, YbarLiqf, rho_liq_own, rho_gas_own);
    cf_nei = thermo.Cfm(p_nei, T_nei, YbarLiqf, rho_liq_nei, rho_gas_nei);

    Info << "max/min cf_own:" 
         << max(cf_own).value() << "/" << min(cf_own).value() << endl;
    Info << "max/min cf_nei:" 
         << max(cf_nei).value() << "/" << min(cf_nei).value() << endl;

    correctCentralACMIInterpolation(cf_nei);    
}


#endif

