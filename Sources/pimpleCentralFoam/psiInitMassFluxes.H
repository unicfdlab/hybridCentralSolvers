W = rho*Rm;

UpdateSonicSpeedOnFaces
(
    thermo,
    own,
    nei,
    p_own,
    p_nei,
    Rmf,
    cf_own,
    cf_nei
);

CalculateKnpWeights
(
    noFluxWeighting,
    fluxScheme,
    v_zero,
    thermo,
    mesh,
    W, // use p/T instead of rho
    U,
    c,
    psi,
    own,
    nei,
    mesh_phi,
    W_own,
    W_nei,
    phiv_own,
    phiv_nei,
    cf_own,
    cf_nei,
    cSf_own,
    cSf_nei,
    ap,
    am,
    alpha_own,
    alpha_nei,
    aSf
);

PsiUpdateCentralFields
(
    noFluxWeighting,
    MRF,
    mesh,
    W, // p/T
    HbyA,
    rAU,
    psi,     // 1/(Rm*T)
    T,
    Rm,
    own,
    nei,
    mesh_phi,
    kappa,
    alpha_own,
    alpha_nei,
    aSf,
    Rmf,
    W_own,
    W_nei,
    phiv_own,
    phiv_nei,
    aphiv_own,
    aphiv_nei,
    phid_own,
    phid_nei,
    psi_own,
    psi_nei,
    psif,
    Dp_own,
    Dp_nei
);