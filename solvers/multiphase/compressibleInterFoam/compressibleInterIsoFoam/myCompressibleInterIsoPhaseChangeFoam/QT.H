//Info<< "Calculating QT" << endl;
volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));
//Info<< "OK limited\n" << endl;

//volScalarField pCoeff(1/rho1 - 1/rho2);

//Info<< "OK pCoeff" << endl;

//volScalarField TSat(1.0/(1.0/T0 - R*log(p/P0)/(W*deltaHvap)));

volScalarField mDotcT
(
    mcCoeff*sqrt(1/TSat)*(1.0-limitedAlpha1)
   *deltaHvap/((1.0/rho2-1.0/rho1)*TSat)
   *mag(fvc::grad(alpha1))*unitLength
   *pos0(T-TSat)
);
volScalarField mDotvT
(
    mvCoeff*sqrt(1/TSat)*limitedAlpha1
   *deltaHvap/((1.0/rho2-1.0/rho1)*TSat)
   *mag(fvc::grad(alpha1))*unitLength 
   *neg(T-TSat)
);

volScalarField QT(deltaHvap*(mDotcT+mDotvT));

Info<< "Heat flux source \n"
    << "Min(QT) = " << min(QT).value()
    << " Max(QT) = " << max(QT).value() << "\n"
    << endl;

/*
volScalarField::Internal vDotcP
(
    IOobject
    (
        "vDotcP",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

volScalarField::Internal vDotvP
(
    IOobject
    (
        "vDotvP",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);


forAll(pCoeff, celli)
{
    vDotvP[celli] = pCoeff[celli]*mcCoeff*limitedAlpha[celli];
    vDotcP[celli] = pCoeff[celli];
}*/

//Info<< "vDotP.H OK\n" << endl;

