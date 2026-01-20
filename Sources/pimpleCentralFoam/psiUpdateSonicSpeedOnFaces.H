#ifndef Psi_Update_Sonic_Speed_On_Faces_H
#define Psi_Update_Sonic_Speed_On_Faces_H

#include "updateSonicSpeedOnFaces.H"

template<>
void UpdateSonicSpeedOnFaces
(
    psiThermo &thermo,
    const surfaceScalarField &own,
    const surfaceScalarField &nei,
    const surfaceScalarField &p_own,
    const surfaceScalarField &p_nei,
    const surfaceScalarField &Rmf,
    surfaceScalarField &cf_own,
    surfaceScalarField &cf_nei
)
{
    const volScalarField& T = thermo.T();

    surfaceScalarField T_own =
        fvc::interpolate(T, own, "reconstruct(psi)");
    surfaceScalarField T_nei =
        fvc::interpolate(T, nei, "reconstruct(psi)");

    const volScalarField Cp (thermo.Cp());
    surfaceScalarField Cp_own =
        fvc::interpolate(Cp, own, "reconstruct(psi)");
    surfaceScalarField Cp_nei =
        fvc::interpolate(Cp, nei, "reconstruct(psi)");
    const volScalarField Cv (thermo.Cv());
    surfaceScalarField Cv_own =
        fvc::interpolate(Cv, own, "reconstruct(psi)");
    surfaceScalarField Cv_nei =
        fvc::interpolate(Cv, nei, "reconstruct(psi)");

    correctCentralACMIInterpolation(T_nei);
    correctCentralACMIInterpolation(Cp_nei);
    correctCentralACMIInterpolation(Cv_nei);

    cf_own = sqrt(Rmf*T_own*Cp_own/Cv_own);
    cf_nei = sqrt(Rmf*T_nei*Cp_nei/Cv_nei);

    correctCentralACMIInterpolation(cf_nei);
}


#endif

