//set correctPhi to true for start time
if (runTime.timeIndex() == runTime.startTimeIndex())
{
    correctPhi = true;
}
else
{
    correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);
}

checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);
