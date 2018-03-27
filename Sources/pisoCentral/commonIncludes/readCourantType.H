word CourantType ("cellCourant");
runTime.controlDict().lookup("CourantType") >> CourantType;
Switch writeCourant(runTime.controlDict().lookupOrDefault<Switch>("writeCourant", false));
