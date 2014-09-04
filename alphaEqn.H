{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phia - phib);

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phic, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, (1.0-alpha), schemer), alpha, schemer)
        );

        alphaEqn.relax();
        alphaEqn.solve();

// in case of exceeding max packing limit
        int i=0;
        while (i<100 && max(alpha).value()>=alphaMax)
        {
             #include "packingLimiter.H"
             i++;
        }
        alpha.min(0.634);
        alpha.max(0);

        beta = scalar(1) - alpha;

        Info<< "Dispersed phase volume fraction = "
            << alpha.weightedAverage(mesh.V()).value()
            << "  Min(alpha) = " << min(alpha).value()
            << "  Max(alpha) = " << max(alpha).value()
            << endl;
    }
}

rho = alpha*rhoa + beta*rhob;
