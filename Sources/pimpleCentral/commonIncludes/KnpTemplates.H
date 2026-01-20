template<class ThermoT>
inline void CalculateKnpWeights
(
    bool noFluxWeighting,
    const Foam::word fluxScheme,
    const dimensionedScalar v_zero,
    const ThermoT& thermo,
    const fvMesh& mesh,
    const volScalarField& W,
    const volVectorField& U,
    volScalarField& c,
    const volScalarField& psi,
    const surfaceScalarField& own,
    const surfaceScalarField& nei,
    const surfaceScalarField& mesh_phi,
    surfaceScalarField& w_own,
    surfaceScalarField& w_nei,
    surfaceScalarField& phiv_own,
    surfaceScalarField& phiv_nei,
    const surfaceScalarField& cf_own,
    const surfaceScalarField& cf_nei,
    surfaceScalarField& cSf_own,
    surfaceScalarField& cSf_nei,
    surfaceScalarField& ap,
    surfaceScalarField& am,
    surfaceScalarField& alpha_own,
    surfaceScalarField& alpha_nei,
    surfaceScalarField& aSf
)
{
    if (noFluxWeighting)
    {
        phiv_own    =
            (fvc::interpolate(U, own, "reconstruct(" + U.name()+")")) & mesh.Sf();
        phiv_nei    =
            (fvc::interpolate(U, nei, "reconstruct(" + U.name()+")")) & mesh.Sf();
    }
    else
    {
        w_own     = fvc::interpolate(W, own, "reconstruct(" + W.name()+")");
        w_nei     = fvc::interpolate(W, nei, "reconstruct(" + W.name()+")");
        correctCentralACMIInterpolation(w_nei);
        phiv_own    =
            (fvc::interpolate(W*U, own, "reconstruct(" + U.name()+")") / w_own) & mesh.Sf();
        phiv_nei    =
            (fvc::interpolate(W*U, nei, "reconstruct(" + U.name()+")") / w_nei) & mesh.Sf();
    }
    
    if (mesh.moving())
    {
        phiv_own   -= mesh_phi;
        phiv_nei   -= mesh_phi;
    }
    
    cSf_own     = cf_own * mesh.magSf();
    cSf_own.setOriented(true);
    cSf_nei     = cf_nei * mesh.magSf();
    cSf_nei.setOriented(true);
    
    ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
    am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
    
    if (fluxScheme == "Kurganov")
    {
        alpha_own   = ap/(ap - am);
        aSf         = am*alpha_own;
        alpha_nei   = 1.0 - alpha_own;
    }
    else if (fluxScheme == "Tadmor")
    {
        aSf         = -0.5*max(mag(ap), mag(am));
    }
    else
    {
        FatalErrorIn
        (
                "updateCentralWeights.H"
        )   << "Unknown central flux scheme " << fluxScheme << nl << nl
            << "Valid flux schemes are:" << nl
            << "1) Kurganov" << nl
            << "2) Tadmor" << endl;
    }
}

template <class MrfT>
inline void PsiUpdateCentralFields
(
    bool noFluxWeighting,
    MrfT& MRF,
    const fvMesh& mesh,
    const volScalarField& W, // p/T
    const volVectorField& HbyA,
    const volScalarField& rAU,
    const volScalarField& psi,     // 1/(Rm*T)
    const volScalarField& T,
    const volScalarField& Rm,
    const surfaceScalarField& own,
    const surfaceScalarField& nei,
    const surfaceScalarField& mesh_phi,
    const surfaceScalarField& kappa,
    const surfaceScalarField& alpha_own,
    const surfaceScalarField& alpha_nei,
    const surfaceScalarField& aSf,
    const surfaceScalarField& Rmf,
    surfaceScalarField& w_own,
    surfaceScalarField& w_nei,
    surfaceScalarField& phiv_own,
    surfaceScalarField& phiv_nei,
    surfaceScalarField& aphiv_own,
    surfaceScalarField& aphiv_nei,
    surfaceScalarField& phid_own,
    surfaceScalarField& phid_nei,
    surfaceScalarField& psi_own,
    surfaceScalarField& psi_nei,
    surfaceScalarField& psif,
    surfaceScalarField& Dp_own,
    surfaceScalarField& Dp_nei
)
{
    psif = linearInterpolate(psi);
    volScalarField invT ("invT", 1.0 / T);
    psi_own     = fvc::interpolate(invT, own, "reconstruct(psi)")/Rmf*kappa
                  + (1.0 - kappa)*psif;
    psi_nei     = fvc::interpolate(invT, nei, "reconstruct(psi)")/Rmf*kappa
                  + (1.0 - kappa)*psif;

    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    if (noFluxWeighting)
    {
        phiv_own    = (fvc::interpolate(HbyA, own, "reconstruct(U)")) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(HbyA, nei, "reconstruct(U)")) & mesh.Sf();
    }
    else
    {
        phiv_own    = (fvc::interpolate(W*HbyA, own, "reconstruct(U)") / w_own) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(W*HbyA, nei, "reconstruct(U)") / w_nei) & mesh.Sf();
    }
    if (mesh.moving())
    {
        phiv_own   -= mesh_phi;
        phiv_nei   -= mesh_phi;
    }
    MRF.makeRelative(phiv_own);
    MRF.makeRelative(phiv_nei);
    aphiv_own   = alpha_own * phiv_own - aSf;
    aphiv_nei   = alpha_nei * phiv_nei + aSf;
    phid_own    = aphiv_own * psi_own;
    phid_nei    = aphiv_nei * psi_nei;
    Dp_own = alpha_own*fvc::interpolate(W * rAU, own, "reconstruct(Dp)")/Rmf;
    Dp_nei = alpha_nei*fvc::interpolate(W * rAU, nei, "reconstruct(Dp)")/Rmf;
}

template <class MrfT>
inline void RhoUpdateCentralFields
(
    bool noFluxWeighting,
    MrfT& MRF,
    const fvMesh& mesh,
    const volScalarField& rho, // rho
    const volScalarField& rhoHat, // rho - psi*p
    const volVectorField& HbyA,
    const volScalarField& rAU,
    const volScalarField& psi,     // 1/(Rm*T)
    const volScalarField& T,
    const surfaceScalarField& own,
    const surfaceScalarField& nei,
    const surfaceScalarField& mesh_phi,
    const surfaceScalarField& kappa,
    const surfaceScalarField& alpha_own,
    const surfaceScalarField& alpha_nei,
    const surfaceScalarField& rho_own,
    const surfaceScalarField& rho_nei,
    const surfaceScalarField& aSf,
    surfaceVectorField& rhoU_own,
    surfaceVectorField& rhoU_nei,
    surfaceScalarField& rhoHatf,
    surfaceScalarField& rhoHat_own,
    surfaceScalarField& rhoHat_nei,
    surfaceScalarField& phiRhoHat_own,
    surfaceScalarField& phiRhoHat_nei,
    surfaceScalarField& phiv_own,
    surfaceScalarField& phiv_nei,
    surfaceScalarField& aphiv_own,
    surfaceScalarField& aphiv_nei,
    surfaceScalarField& phid_own,
    surfaceScalarField& phid_nei,
    surfaceScalarField& psi_own,
    surfaceScalarField& psi_nei,
    surfaceScalarField& psif,
    surfaceScalarField& Dp_own,
    surfaceScalarField& Dp_nei
)
{
    psif = linearInterpolate(psi);
    psi_own       = fvc::interpolate(psi, own, "reconstruct(psi)")*kappa
                      + (1.0 - kappa)*psif;
    psi_nei       = fvc::interpolate(psi, nei, "reconstruct(psi)")*kappa
                      + (1.0 - kappa)*psif;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    
    rhoHatf = linearInterpolate(rhoHat);
    rhoHat_own    = fvc::interpolate(rhoHat, own, "reconstruct(rhoHat)")*kappa
                      + (1.0 - kappa)*rhoHatf;
    rhoHat_nei    = fvc::interpolate(rhoHat, nei, "reconstruct(rhoHat)")*kappa
                      + (1.0 - kappa)*rhoHatf;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(rhoHat_nei);
    
    if (noFluxWeighting)
    {
        phiv_own    = (fvc::interpolate(HbyA, own, "reconstruct(U)")) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(HbyA, nei, "reconstruct(U)")) & mesh.Sf();
    }
    else
    {
        rhoU_own      = fvc::interpolate(rho*HbyA, own, "reconstruct(U)");
        rhoU_nei      = fvc::interpolate(rho*HbyA, nei, "reconstruct(U)");
        phiv_own      = (rhoU_own / rho_own) & mesh.Sf();
        phiv_nei      = (rhoU_nei / rho_nei) & mesh.Sf();
    }
    
    if (mesh.moving())
    {
        phiv_own   -= mesh_phi;
        phiv_nei   -= mesh_phi;
    }
    
    MRF.makeRelative(phiv_own);
    MRF.makeRelative(phiv_nei);
    
    aphiv_own     = alpha_own * phiv_own - aSf;
    aphiv_nei     = alpha_nei * phiv_nei + aSf;
    
    phid_own      = aphiv_own * psi_own;
    phid_nei      = aphiv_nei * psi_nei;
    
    phiRhoHat_own = aphiv_own * rhoHat_own;
    phiRhoHat_nei = aphiv_nei * rhoHat_nei;
    
    Dp_own = alpha_own*fvc::interpolate(rho * rAU, own, "reconstruct(Dp)");
    Dp_nei = alpha_nei*fvc::interpolate(rho * rAU, nei, "reconstruct(Dp)");
}

inline void BlendWithLinear
(
    const surfaceScalarField &kappa,
    const surfaceScalarField &onemkappa,
    surfaceScalarField &xf,
    surfaceScalarField &x_own,
    surfaceScalarField &x_nei
) {
    xf    *= onemkappa,
    x_own *= kappa;
    x_own += xf;
    x_nei *= kappa;
}
