scalar pi = constant::mathematical::pi;

// following is general form of mean flow + sine wave and 2nd order stokes waves
gradPOSC = gradPMEAN+gradPAMP1*Foam::cos(initTheta*(pi)/180.0 + 2.0*(pi)*runTime.value()/oscpT)+gradPAMP2*Foam::sin(initTheta*(pi)/90.0 + 4.0*(pi)*runTime.value()/oscpT);     // Add oscilitary pressure gradient 

fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);

{
    {
        rho = alpha*rhoa+beta*rhob;
// defined in callKineticTheory.H
        gradUaT = fvc::grad(Ua)().T();
        
        volTensorField Rca
        (
            "Rca",
            ((2.0/3.0)*I)*(sqr(Ct)*k + nuEffa*tr(gradUaT)) - nuEffa*gradUaT
        );

        if (kineticTheory.on())
        {
            Rca -= ((lambdaUa/((alpha+1e-8)*rhoa))*tr(gradUaT))*tensor(I);
        }

        surfaceScalarField phiRa =
            -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
            /fvc::interpolate((alpha+ scalar(1e-8)));
        volVectorField susa
        (
            "susa",
// turbulent suspension in all directions
            SUS*(K/rhoa)*nutb*(beta/((alpha+scalar(1e-08)))*fvc::grad(alpha))
        );

        UaEqn =
        (
            (scalar(1.0) + Cvm*rhob*beta/rhoa)*
            (
                fvm::ddt(Ua)
              + fvm::div(phia, Ua, "div(phia,Ua)")
              - fvm::Sp(fvc::div(phia), Ua)
            )
         
// effective stress term from dispersed phase
          - fvm::laplacian(nuEffa, Ua)
          + fvc::div(Rca)
          + fvm::div(phiRa, Ua, "div(phiRa,Ua)")
          - fvm::Sp(fvc::div(phiRa), Ua)
          + (fvc::grad(alpha)/((alpha+ scalar(1e-08))) & Rca)
         ==
//        - gradPs                   //the particle granular pressure is moved to p-equation
//        - gradPf                   //the particle frictional pressure is moved to p-equation
          - fvm::Sp(beta*K/rhoa, Ua) // treat drag term implicitly
//        + beta/rhoa*K*Ub            // Explicit drag transfered to p-equation
          - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)
//        - susa
//        + gradPOSC/(tilt*rhob+(1.0-tilt)*rhoa) // moved to p-equation Add oscilitary pressure gradient
//        + g                      // Buoyancy term transfered to p-equation
        );
          
        UaEqn.relax();
    }

    {
        volTensorField gradUbT = fvc::grad(Ub)().T();
        volTensorField Rcb
        (
            "Rcb",
// move k to p equation, to avoid BC oscillation
            ((2.0/3.0)*I)*(nuEffb*tr(gradUbT)) - nuEffb*gradUbT
        );

        surfaceScalarField phiRb =
            -fvc::interpolate(nuEffb)*mesh.magSf()*(fvc::snGrad(beta))
            /fvc::interpolate(beta);

        volVectorField susb
        (
            "susb",
            SUS*(K/rhob)*nutb*(fvc::grad(alpha))
        );

       UbEqn =
        (
            (scalar(1.0) + Cvm*rhob*alpha/rhob)*
            (
                fvm::ddt(Ub)
              + fvm::div(phib, Ub, "div(phib,Ub)")
              - fvm::Sp(fvc::div(phib), Ub)
            )

          - fvm::laplacian(nuEffb, Ub)
          + fvc::div(Rcb)
          + fvm::div(phiRb, Ub, "div(phiRb,Ub)")
          - fvm::Sp(fvc::div(phiRb), Ub)
          + (fvc::grad(beta)/beta & Rcb)
         ==
          - fvm::Sp(alpha*K/rhob, Ub)   // treat drag term implicitly 
//          + alpha/rhob*K*Ua              // Explicit drag transfered to p-equation
          + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
//        + susb                          // add the suspension to p-equation
//        + (scalar(1.0)/rhob)*gradPOSC  // driving force moved to p-equation
//        + g                       // Buoyancy term transfered to p-equation
        );
        UbEqn.relax();
    }
}
