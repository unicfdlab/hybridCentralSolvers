{
    aphiv_own = phi_own / rho_own;
    aphiv_nei = phi_nei / rho_nei;
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    amaxSf = max(mag(aphiv_own), mag(aphiv_nei)) + amaxSmall;
    amaxSf.setOriented(true);
    
    kappaFuncPtr().update();
    kappa = kappaFuncPtr().kappa()();

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phi_own += (1.0 - kappa) * phi_nei;
    phi_nei = kappa * phi_nei;
}

