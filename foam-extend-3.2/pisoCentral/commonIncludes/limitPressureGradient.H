{
    limGradP = fvc::grad(p);
    forAll(badQualityCells, iCell)
    {
	gradp[badQualityCells[iCell]] = limGradP[badQualityCells[iCell]];
    }
}

