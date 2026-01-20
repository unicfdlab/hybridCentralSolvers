List<scalar> ScNumbers (composition.species().size(), 1.0);

if (thermo.found("ScNumbers"))
{
    List<Tuple2<word,scalar> > ScNumbersTuples;
    thermo.lookup("ScNumbers") >> ScNumbersTuples;
    
    forAll(ScNumbersTuples, iTuple)
    {
        word cmptName = ScNumbersTuples[iTuple].first();
        scalar cmptSc   = ScNumbersTuples[iTuple].second();
        if (composition.species().found(cmptName))
        {
            label cmptId = composition.species()[cmptName];
            ScNumbers[cmptId] = cmptSc;
        }
    }
}
