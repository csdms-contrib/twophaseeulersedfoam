if(turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }
// to be consistent, let's define nutb = C_mu*k_f^2/epsilon_f
    tmp<volTensorField> tgradUb = fvc::grad(Ub);
    volTensorField TgradUb = fvc::grad(Ub)().T();
// same as that in momentum equation
    volScalarField G = fvc::grad(Ub) && (-(2.0/3.0)*(k + nutb*tr(TgradUb))*I + nutb*(tgradUb+TgradUb));
    tgradUb.clear();
    dimensionedScalar Tpsmall
       (
           "Tpsmall",
           dimensionSet(1, -3, -1, 0, 0, 0, 0),
           scalar(1e-40)
        );                                  
    dimensionedScalar ksmall
    (
           "ksmall",
           dimensionSet(0, 2, -2, 0, 0, 0, 0),
           scalar(1e-20)
    );
    dimensionedScalar epsilonsmall
    (
           "epsilonsmall",
           dimensionSet(0, 2, -3, 0, 0, 0, 0),
           scalar(1e-20)
    );
// found that drag coefficient beta = K, see T.J.Hsu 2004, fixed the definition of K already,namely, rhos/Tp=(1-alpha)*K
//    Tp=rhoa/(beta*K+Tpsmall); 
//    Tl = k/(6.0*epsilon);
//    tmf = Foam::exp(-B*Tp/Tl);
    tmf = Foam::exp(-B*rhoa*6.0*epsilon/(k*max(beta*K,Tpsmall)));
// just need the horizontal ones;
    volScalarField rtemp = (Ub-Ua) & g;
    volVectorField Urtemp = rtemp*g/magSqr(g);
    volVectorField Urh = (Ub-Ua) - Urtemp;
// ESD1 for horizontal, ESD4 for vertical density stratification    
// for 1D ESD1 ==0, obviously
    volScalarField ESD1 = nutb*(K/rhob)*(fvc::grad(alpha)) & Urh;
    volScalarField ESD2 = -2.0*K*alpha*(1.0-tmf)/rhob; //k_f was included in the equation implicitly
    volScalarField ESD3 = K*alpha*magSqr(Ub-Ua)/(beta*rhob);
    volScalarField ESD4 = -(nutb/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha);
    ESD4 = ESD4 + KE1*ESD1;
// the following term is implemented in epsilon equation.
    volScalarField ESD5 = -(nutb/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha);
// if ESD5 > 0, then it should be implemented, else, no
    ESD5 = max(ESD5,scalar(0.0)*ESD5);
    ESD5 = ESD5 + KE1*ESD1;
    volScalarField ESD = KE3*ESD3;   //need to check here!
    
    #include "wallFunctions.H"

    // Dissipation equation
    //epsilon = epsilon_f
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon)
      + fvm::div(phib, epsilon)
      - fvm::Sp(fvc::div(phib),epsilon)
      - fvm::laplacian
        (
            (alphaEps*nutb+nub), epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
      ==
         fvm::SuSp(C1*G/k ,epsilon)
//         C1*G*epsilon/k
// The following one line should be Sp for stability
       + fvm::Sp(-C2*epsilon/k, epsilon)
       + fvm::Sp(C3ep*(ESD/k),epsilon)
       + fvm::Sp(C3ep*KE2*ESD2,epsilon)
       + fvm::Sp((C4ep*KE4*ESD5/k),epsilon)
    );

    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-20));

    // Turbulent kinetic energy equation
    // k=k_f 
    fvScalarMatrix kEqn
    (
        fvm::ddt(k)
      + fvm::div(phib, k)
      - fvm::Sp(fvc::div(phib),k)
      - fvm::laplacian
        (
            (alphak*nutb+nub), k,
            "laplacian(DkEff,k)"
        )
      ==
        fvm::SuSp(G/k,k)
//        G
// The following one line should be Sp for stability
      + fvm::Sp(-epsilon/k,k)
      + fvm::Sp(ESD/k,k)
      + fvm::Sp(KE4*ESD4/k,k)
      + fvm::Sp(KE2*ESD2,k)
    );
    kEqn.relax();
    kEqn.solve();

    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-20));

    //- Re-calculate turbulence viscosity
    nutb = Cmu*sqr(k)/epsilon;
 
    #include "wallViscosity.H"
}
Info<< "max(nutb) = " << max(nutb).value() <<" min(nutb) = "<<min(nutb).value()<< endl;
// limit nutb for stability reason, 
nutb.min(dimensionedScalar("nutlimit", nub.dimensions(), 1));
nuEffb = (nutb + nub);
