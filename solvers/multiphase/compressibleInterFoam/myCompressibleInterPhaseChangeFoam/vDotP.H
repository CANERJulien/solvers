//Info<< "Calculating vDotP" << endl;
volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));
//Info<< "OK limited\n" << endl;

/*
Pair<tmp<volScalarField>> mDotP
                            (
                                mcCoeff*(1.0-limitedAlpha1)*pos0(p - pSat),
                                (-mvCoeff)*limitedAlpha1*neg(p - pSat)
                            );
Info<< "OK mDotP\n" << endl;

//dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
const volScalarField& pCoeff(1.0/rho1 - 1.0/rho2);
Info<< "OK pCoeff" << endl;
Pair<tmp<volScalarField>> vDotP (pCoeff*mDotP[0], pCoeff*mDotP[1]);
Info<< "OK vDotP" << endl;
const volScalarField& vDotcP = vDotP[0];
const volScalarField& vDotvP = vDotP[1];
*/

volScalarField pCoeff(1/rho1 - 1/rho2);

//Info<< "OK pCoeff" << endl;

//volScalarField vDotcP(pCoeff*mcCoeff*(1.0-limitedAlpha1)*pos0(p-pSat));
volScalarField vDotcP
(
    pCoeff*mcCoeff*sqrt(1/T)*(1.0-limitedAlpha1)
//   *mag(fvc::grad(alpha1))*unitLength
   *pos0(p-pSat)
);
//volScalarField vDotvP(pCoeff*(-mvCoeff)*limitedAlpha1*neg(p-pSat));
volScalarField vDotvP
(
    pCoeff*(-mvCoeff*sqrt(1/T))*limitedAlpha1
//   *mag(fvc::grad(alpha1))*unitLength
   *neg(p-pSat)
);

Info<< "Test Volumetric condensation and vaporization rates \n"
    << "Min(vDotcP) = " << min(vDotcP).value()
    << " Max(vDotcP) = " << max(vDotcP).value() << "\n"
    << "Min(vDotvP) = " << min(vDotvP).value()
    << " Max(vDotvP) = " << max(vDotvP).value()
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

// Ecrire le debit massique
/*
if (mesh.time().outputTime())
{
    volScalarField vDotv
    (
    "vDotv",
    vDotvP*(p-pSat)
    );
    volScalarField vDotc
    (
    "vDotc",
    vDotcP*(p-pSat)
    );



    vDotv.write();
    vDotc.write();
}
*/
