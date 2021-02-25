k=k1*min(alpha1,scalar(1))+k2*max(alpha2,scalar(0));
rhocp=(rho1cp1)*alpha1+(rho2cp2)*alpha2;       
//Q=m*lv*mag(fvc::grad(alpha1))*unitLength;


Pair<tmp<volScalarField>> mDotAlphal =
        mixture->mDotAlphal();

//volScalarField Q=(mDotAlphal[0]+mDotAlphal[1])*lv*mag(fvc::grad(alpha1))*unitLength;
Q=(mDotAlphal[0]+mDotAlphal[1])*lv*mag(fvc::grad(alpha1))*unitLength;

Info<< "Ici 40 = " << runTime.timeName() << nl << endl;

surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

surfaceScalarField phir("phir", phic*interface.nHatf());

//surfaceScalarField phirAlpha(phir*fvc::interpolate(alpha2));

phirAlpha=phir*fvc::interpolate(alpha2);

//phi1=Kc*interface.nHatf()*max(phi*interface.nHatf()/sqr(mesh.magSf()));
Info<< "Ici 41 = " << runTime.timeName() << nl << endl;
//word alpharScheme("div(phirb,alpha)");

fvScalarMatrix Teqn
 (
   (rhocp)*(fvm::ddt(T)
            +fvm::div(phi,T)
            +fvm::div(phirAlpha,T)
           ) 
               

    // +fvm::div(phir,(1-alpha1)*T)
    -fvm::laplacian(k,T)
    -Q
    
  );
 Info<< "Ici 42 = " << runTime.timeName() << nl << endl;
  Teqn.solve();
