// This version is modified so that Ua=0 in flow direction when alpha>alphaMax

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

    forAll (alpha, celll)
    {
       const cell & cc = mesh.cells()[celll];
       labelList pLabels(cc.labels(ff));
       pointField pLocal(pLabels.size(), vector::zero);

       forAll (pLabels, pointi)
       pLocal[pointi] = pp[pLabels[pointi]];

       //scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));

       if(alpha[celll]>=alphaMinFriction.value())
       {
// in this way, the fluid veocity at the bottom would not be affected by the wall function.
           Ua[celll] = 0*Ua[celll];
           Ub[celll] = 0*Ub[celll];
       }
    }
phia = fvc::interpolate(Ua) & mesh.Sf();
phib = fvc::interpolate(Ub) & mesh.Sf();
phi = alphaf*phia + betaf*phib;
