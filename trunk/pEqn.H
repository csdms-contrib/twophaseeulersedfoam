{
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField betaf(scalar(1) - alphaf);
 
    volScalarField rUaA(1.0/UaEqn.A());
    volScalarField rUbA(1.0/UbEqn.A());

    phia == (fvc::interpolate(Ua) & mesh.Sf());
    phib == (fvc::interpolate(Ub) & mesh.Sf());

    rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf(fvc::interpolate(rUbA));

    surfaceScalarField phiDraga
    (
       fvc::interpolate(beta/rhoa*K*rUaA)*phib -
       (fvc::interpolate(SUS*(K/rhoa)*nutb*(beta/(alpha+scalar(1e-08))))*fvc::snGrad(alpha))*rUaAf*mesh.magSf() -
       (fvc::interpolate(scalar(1.0)/(rhoa*((alpha+scalar(1e-08)))))*fvc::snGrad(kineticTheory.pa()))*rUaAf*mesh.magSf() -
       (fvc::interpolate(scalar(1.0)/(rhoa*((alpha+scalar(1e-08)))))*fvc::snGrad(kineticTheory.pf_()))*rUaAf*mesh.magSf()+
       rUaAf*(g & mesh.Sf())+
       rUaAf*((gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) & mesh.Sf())
    );
       
    surfaceScalarField phiDragb
    (
       fvc::interpolate(alpha/rhob*K*rUbA)*phia +
       (fvc::interpolate(SUS*(K/rhob)*nutb)*fvc::snGrad(alpha))*rUbAf*mesh.magSf() -
       fvc::interpolate(scalar(2.0)/(beta*scalar(3.0)))*fvc::snGrad(k*beta)*rUbAf*mesh.magSf() +
       rUbAf*(g & mesh.Sf())+
       rUbAf*((gradPOSC/rhob) & mesh.Sf())
    );

    Ua = rUaA*UaEqn.H();
    Ub = rUbA*UbEqn.H();

    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (isA<symmetryFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (p.boundaryField().types()[patchi] == "symmetryPlane")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }

    surfaceScalarField phiUa = (fvc::interpolate(Ua) & mesh.Sf()) + fvc::ddtPhiCorr(rUaA, Ua, phia);
    surfaceScalarField phiUb = (fvc::interpolate(Ub) & mesh.Sf()) + fvc::ddtPhiCorr(rUbA, Ub, phib);

    adjustPhi(phiUa, Ua, p);
    adjustPhi(phiUb, Ub, p);

    phia = phiUa + phiDraga;
    phib = phiUb + phiDragb;

    phi = alphaf*phia + betaf*phib;
 
    surfaceScalarField Dp
    (
        "(rho*(1|A(U)))",
        alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        p.correctBoundaryConditions();

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phia -= rUaAf*SfGradp/rhoa;
            phib -= rUbAf*SfGradp/rhob;
            phi -= pEqn.flux(); 
        }
    }
    Ua += rUaA*fvc::reconstruct((phia - phiUa)/rUaAf);
    Ua.correctBoundaryConditions();
    Ub += rUbA*fvc::reconstruct((phib - phiUb)/rUbAf);
    Ub.correctBoundaryConditions();
    U = alpha*Ua + beta*Ub;
    #include "continuityErrs.H"
}
