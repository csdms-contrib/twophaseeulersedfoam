Info<< "\nReading transportProperties for gradP\n" <<endl;

    dimensionedScalar tilt
    (
        transportProperties.lookup("tilt")
    );

    dimensionedVector gradPAMP1
    (
        transportProperties.lookup("gradPAMP1")
    );
    dimensionedVector gradPAMP2
    (
        transportProperties.lookup("gradPAMP2")
    );

    dimensionedVector gradPMEAN
    (
        transportProperties.lookup("gradPMEAN")
    );

    dimensionedVector gradPOSC
    (
        transportProperties.lookup("gradPOSC")
    );

      dimensionedScalar oscpT
    (
        transportProperties.lookup("oscpT")
    );
      dimensionedScalar initTheta
    (
        transportProperties.lookup("initTheta")
    );
      dimensionedScalar SUS
    (
        transportProperties.lookup("SUS")
    );
      dimensionedScalar C3ep
    (
        transportProperties.lookup("C3ep")
    );
      dimensionedScalar C4ep
    (
        transportProperties.lookup("C4ep")
    );
            dimensionedScalar KE1
    (
        transportProperties.lookup("KE1")
    );  
            dimensionedScalar KE2
    (
        transportProperties.lookup("KE2")
    );
            dimensionedScalar KE3
    (
        transportProperties.lookup("KE3")
    );
    dimensionedScalar KE4
    (
        transportProperties.lookup("KE4")
    );
    dimensionedScalar B
    (
        transportProperties.lookup("B")
    );
    dimensionedScalar alphaMinFriction
    (
        transportProperties.lookup("alphaMinFriction")
    );
      Info<< "Initializing with specified pressure gradient:"<<gradPOSC.value()<<endl
            << endl;
