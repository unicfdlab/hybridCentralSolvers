{
    Info << "Marking cells with bad quality to limit pressure gradient" << endl;
    
    badQualityCells.resize(0);
    
    const scalar nonOrthogonalityThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("nonOrthogonalityThreshold")
        )
    );
    
    const scalar skewnessThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("skewnessThreshold")
        )
    );

    cellQuality cq (mesh);
    
    scalarField cellNonOrtho (cq.nonOrthogonality());
    
    scalarField cellSkewness (cq.skewness());
    
    forAll(cellNonOrtho, iCell)
    {
        if ( (cellNonOrtho[iCell] >= nonOrthogonalityThreshold) ||
            (cellSkewness[iCell] >= skewnessThreshold) )
        {
            badQualityCells.append(iCell);
            
        }
    }
}
