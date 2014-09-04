    if (runTime.outputTime())
    {
        volVectorField Ur
        (
            IOobject
            (
                "Ur",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ub - Ua
        );
        volVectorField susb 
        (
            IOobject
            (
                "susb",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            SUS*((draga->K(mag(Ur)))/rhob)*nutb*(fvc::grad(alpha))
//            SUS*((draga->K(mag(Ur)))/rhob)*nutb*(fvc::grad(alpha) & g)*g/magSqr(g)
        );

        volTensorField Tauf
        (
           IOobject
           (
                "Tauf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
           ),
           nuEffb*(fvc::grad(Ub)+fvc::grad(Ub)().T())
        );
        volTensorField Taus
        (
           IOobject
           (
                "Taus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
           ),
           nuEffa*(fvc::grad(Ua)+fvc::grad(Ua)().T())
        );

        runTime.write();
    }
