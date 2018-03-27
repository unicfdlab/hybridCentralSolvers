/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURownE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

CoNum = 0.0;
meanCoNum = 0.0;
velMag = 0.0;
{
    autoPtr<volScalarField> localCoPtr;
    
    if ( CourantType == "cellCourant" )
    {
        if (mesh.nInternalFaces())
        {
            localCoPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Co",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar("localCo", dimless,  VSMALL),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            
            localCoPtr().primitiveFieldRef() = 
            (
                fvc::surfaceSum(mag(phi))().primitiveField()
                /
                rho.primitiveField()
                /
                mesh.V().field()
            ) * 0.5 * runTime.deltaT().value();
            
            CoNum = gMax(localCoPtr().internalField());
        }
        
        Info<< "Courant Number max: " <<  CoNum << endl;
    }
    else if ( CourantType == "faceCourant" )
    {
        if (mesh.nInternalFaces())
        {
        
            if (writeCourant)
            {
                localCoPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "Co",
                            runTime.timeName(),
                            mesh
                        ),
                        mesh,
                        dimensionedScalar("localCo", dimless,  VSMALL),
                        calculatedFvPatchScalarField::typeName
                    )
                );
            }
            
            surfaceScalarField magUf =
            mag
            (
                phi
                / 
                (
                    rho_own * alpha_own
                    +
                    rho_nei * alpha_nei
                ) / tuMagSf()
            );
            
            surfaceScalarField faceCo =
                mesh.surfaceInterpolation::deltaCoeffs()
                *
                magUf
                *
                runTime.deltaT() * 2.0;
            
            //update boundary Courant numbers
            forAll (faceCo.boundaryField(), iPatch)
            {
                forAll(faceCo.boundaryField()[iPatch], iFace)
                {
                    faceCo.boundaryFieldRef()[iPatch][iFace] /= 2.0;
                }
            }
        
            CoNum = max(faceCo).value();
        
            velMag = max(magUf).value();
            
            if (writeCourant)
            {
                forAll(mesh.cells(), iCell)
                {
                    scalar maxCoCell =  -0.01;
                    
                    const labelList& cellFaces = mesh.cells()[iCell];
                    
                    forAll(cellFaces, iFace)
                    {
                        if (mesh.isInternalFace(cellFaces[iFace]))
                        {
                            if (faceCo[cellFaces[iFace]] > maxCoCell)
                            {
                                maxCoCell = faceCo[cellFaces[iFace]];
                            }
                       }
                    }
                    
                    localCoPtr()[iCell] = maxCoCell;
                }
                
                forAll(localCoPtr().boundaryField(), iPatch)
                {
                    forAll(localCoPtr().boundaryField()[iPatch], iFace)
                    {
                        localCoPtr().boundaryFieldRef()[iPatch][iFace] = 
                            faceCo.boundaryField()[iPatch][iFace];
                    }
                }
            } //if writeCourant
        } //if mesh.nInternalFaces() > 9
    
        Info<< "Courant Number max : " << CoNum
            << " velocity magnitude: " << velMag
            << endl;
    }
    /*
    else if ( CourantType == "mixedFaceCourant" )
    {

        if (mesh.nInternalFaces())
        {
        
            if (writeCourant)
            {
                localCoPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "Co",
                            runTime.timeName(),
                            mesh
                        ),
                        mesh,
                        dimensionedScalar("localCo", dimless,  1.0)
                    )
                );
            }
            
            surfaceScalarField magUf =
            kappa*amaxSf/tuMagSf()
            +
            (1.0 - kappa)*
            mag
            (
                phi
                / 
                (
                    rho_own * alpha_own
                    +
                    rho_nei * alpha_nei
                ) / tuMagSf()
            );
            
            surfaceScalarField faceCo =
                mesh.surfaceInterpolation::deltaCoeffs()
                *
                magUf
                *
                runTime.deltaT() * 2.0;
        
            CoNum = max(faceCo).value();
        
            velMag = max(magUf).value();
            
            if (writeCourant)
            {
                forAll(mesh.cells(), iCell)
                {
                    scalar maxCoCell =  -0.01;
                    
                    const labelList& cellFaces = mesh.cells()[iCell];
                    
                    forAll(cellFaces, iFace)
                    {
                        if (mesh.isInternalFace(cellFaces[iFace]))
                        {
                            if (faceCo[cellFaces[iFace]] > maxCoCell)
                            {
                                maxCoCell = faceCo[cellFaces[iFace]];
                            }
                       }
                    }
                    
                    localCoPtr()[iCell] = maxCoCell;
                }
                
                forAll(localCoPtr().boundaryField(), iPatch)
                {
                    forAll(localCoPtr().boundaryField()[iPatch], iFace)
                    {
                        localCoPtr().boundaryFieldRef()[iPatch][iFace] = 
                            faceCo.boundaryField()[iPatch][iFace];
                    }
                }
            } //if writeCourant
        } //if mesh.nInternalFaces() > 9
    
        Info<< "Courant Number max : " << CoNum
            << " velocity magnitude: " << velMag
            << endl;
    }
    */
    else
    {
        FatalErrorIn
        (
            "compressibleCourantNo.H:"
        )   << "Wrong type of Courant criterion: " << CourantType
        << endl << " must be one of:" 
        << endl << "1) cellCourant"
        << endl << "2) faceCourant"
        //<< endl << "3) mixedFaceCourant"
        << endl << abort(FatalError);
    }
    
    if (localCoPtr.valid() && writeCourant)
    {
        if (runTime.outputTime())
        {
            localCoPtr().write();
        }
    }
}

// ************************************************************************* //
