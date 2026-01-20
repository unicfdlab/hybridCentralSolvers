#ifndef Update_Sonic_Speed_On_Faces_H
#define Update_Sonic_Speed_On_Faces_H

template<typename ThermoT>
void UpdateSonicSpeedOnFaces
(
    ThermoT &thermo,
    const surfaceScalarField &own,
    const surfaceScalarField &nei,
    const surfaceScalarField &p_own,
    const surfaceScalarField &p_nei,
    const surfaceScalarField &Rmf,
    surfaceScalarField &cf_own,
    surfaceScalarField &cf_nei
) {
    volScalarField c = sqrt(thermo.Cp()/thermo.Cv() / thermo.psi());
    cf_own           = fvc::interpolate(c, own, "reconstruct(psi)");
    cf_nei           = fvc::interpolate(c, nei, "reconstruct(psi)");
    Info << "max/min c: " << max(c).value() << "/" << min(c).value()
         << ", default interpolation" << endl;
}
// {
//     const volScalarField& T = thermo.T();
//     basicSpecieMixture& composition = thermo.composition();

//     surfaceScalarField T_own =
//         fvc::interpolate(T, own, "reconstruct(psi)");
//     surfaceScalarField T_nei =
//         fvc::interpolate(T, nei, "reconstruct(psi)");

//     scalar Cpm_own = 0.0;
//     scalar Cvm_own = 0.0;
//     scalar Cpm_nei = 0.0;
//     scalar Cvm_nei = 0.0;
//     PtrList<volScalarField>& Y = composition.Y();
//     PtrList<surfaceScalarField> Yf(Y.size());
//     forAll(Y, ispecie)
//     {
//         Yf.set
//         (
//             ispecie,
//             linearInterpolate(Y[ispecie])
//         ); // TVD-scheme ?
//     }

//     forAll(p_own.primitiveField(), iface)
//     {
//         Cpm_own = 0.0;
//         Cvm_own = 0.0;
//         Cpm_nei = 0.0;
//         Cvm_nei = 0.0;
//         forAll(Yf, ispecie)
//         {
//             Cpm_own += thermo.composition().Cp
//             (
//                 ispecie,
//                 p_own.primitiveField()[iface],
//                 T_own.primitiveField()[iface]
//             )*Yf[ispecie][iface];
//             Cvm_own += thermo.composition().Cv
//             (
//                 ispecie,
//                 p_own.primitiveField()[iface],
//                 T_own.primitiveField()[iface]
//             )*Yf[ispecie][iface];
//             Cpm_nei += thermo.composition().Cp
//             (
//                 ispecie,
//                 p_nei.primitiveField()[iface],
//                 T_nei.primitiveField()[iface]
//             )*Yf[ispecie][iface];
//             Cvm_nei += thermo.composition().Cv
//             (
//                 ispecie,
//                 p_nei.primitiveField()[iface],
//                 T_nei.primitiveField()[iface]
//             )*Yf[ispecie][iface];
//         }
//         cf_own.primitiveFieldRef()[iface] =
//             Foam::sqrt
//             (
//                 (Cpm_own/Cvm_own)*
//                 Rmf.primitiveField()[iface]*
//                 T_own.primitiveField()[iface]
//             );
//         cf_nei.primitiveFieldRef()[iface] =
//             Foam::sqrt
//             (
//                 (Cpm_nei/Cvm_nei)*
//                 Rmf.primitiveField()[iface]*
//                 T_nei.primitiveField()[iface]
//             );
//     }

//     forAll(p_own.boundaryField(), ipatch)
//     {
//         const fvsPatchField<scalar> &pp_own = p_own.boundaryField()[ipatch];
//         const fvsPatchField<scalar> &pT_own = T_own.boundaryField()[ipatch];
//         const fvsPatchField<scalar> &pp_nei = p_nei.boundaryField()[ipatch];
//         const fvsPatchField<scalar> &pT_nei = T_nei.boundaryField()[ipatch];
//         fvsPatchField<scalar> &pcf_own = cf_own.boundaryFieldRef()[ipatch];
//         fvsPatchField<scalar> &pcf_nei = cf_nei.boundaryFieldRef()[ipatch];
//         const fvsPatchField<scalar> &pRmf = Rmf.boundaryField()[ipatch];
//         forAll(pp_own, iface)
//         {
//             Cpm_own = 0.0;
//             Cvm_own = 0.0;
//             Cpm_nei = 0.0;
//             Cvm_nei = 0.0;
//             forAll(Yf, ispecie)
//             {
//                 Cpm_own += thermo.composition().Cp
//                 (
//                     ispecie, pp_own[iface], pT_own[iface]
//                 )*Yf[ispecie].boundaryField()[ipatch][iface];
//                 Cvm_own += thermo.composition().Cv
//                 (
//                     ispecie, pp_own[iface], pT_own[iface]
//                 )*Yf[ispecie].boundaryField()[ipatch][iface]; 
//                 Cpm_nei += thermo.composition().Cp
//                 (
//                     ispecie, pp_nei[iface], pT_nei[iface]
//                 )*Yf[ispecie].boundaryField()[ipatch][iface];
//                 Cvm_nei += thermo.composition().Cv
//                 (
//                     ispecie, pp_nei[iface], pT_nei[iface]
//                 )*Yf[ispecie].boundaryField()[ipatch][iface];
//             }
//             pcf_own[iface] =
//                 Foam::sqrt((Cpm_own/Cvm_own)*pRmf[iface] * pT_own[iface]);
//             pcf_nei[iface] =
//                 Foam::sqrt((Cpm_nei/Cvm_nei)*pRmf[iface] * pT_nei[iface]);
//         }
//     }
//     correctCentralACMIInterpolation(cf_nei);
// }

#endif

