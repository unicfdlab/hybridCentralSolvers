#ifndef correctCentralACMIInterpolation_H
#define correctCentralACMIInterpolation_H

#include "surfaceFieldsFwd.H"

namespace Foam
{

void correctCentralACMIInterpolation(surfaceScalarField& nei_field);

void correctCentralACMIInterpolation(surfaceVectorField& nei_field);

}

#endif

//
//END-OF-FILE
//

