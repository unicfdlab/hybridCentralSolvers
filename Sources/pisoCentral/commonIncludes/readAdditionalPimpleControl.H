{
    pimple.dict().lookup("updateEnergyInPISO") >> updateEnergyInPISO;
    
    pimple.dict().lookup("fluxScheme") >> fluxScheme;

    if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
    {
    }
    else
    {
        FatalErrorIn
        (
        "readAdditionalPimpleControl.H "
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov"
            << abort(FatalError);
    }
    
    if (pimple.dict().found("noFluxWeighting"))
    {
        pimple.dict().lookup("noFluxWeighting") >> noFluxWeighting;
    }
}
