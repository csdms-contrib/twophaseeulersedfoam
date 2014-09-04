volTensorField gradUaT = fvc::grad(Ua)().T();
if (kineticTheory.on())
{
   kineticTheory.solve(gradUaT,k,epsilon,nutb,B,runTime);
   nuEffa = kineticTheory.mua()/((alpha+1e-5)*rhoa);
   lambdaUa = kineticTheory.lambda();
}
   else               // If not using kinetic theory is using Ct model
{
   nuEffa = sqr(Ct)*nutb + nua;
}
if (kineticTheory.on())
{
   gradPs = (fvc::grad(kineticTheory.pa()))/(rhoa*((alpha+scalar(1e-05))));
   gradPf = (fvc::grad(kineticTheory.pf_()))/(rhoa*((alpha+scalar(1e-05))));
   Info<< "Max gradPs =" << max(mag(gradPs)).value() <<",Min gradPs =" << min(mag(gradPs)).value()<<endl;
   Info<< "Max gradPf =" << max(mag(gradPf)).value() <<",Min gradPf =" << min(mag(gradPf)).value()<<endl;
}
