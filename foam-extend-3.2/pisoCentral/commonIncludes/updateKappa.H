{
    aphiv_pos = phi_pos / rho_pos;
    aphiv_neg = phi_neg / rho_neg;
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg)) + amaxSmall;
    
    kappaFuncPtr().update();
    kappa = kappaFuncPtr().kappa()();

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phi_pos = phi_pos + (1.0 - kappa) * phi_neg;
    phi_neg = kappa * phi_neg;
}

