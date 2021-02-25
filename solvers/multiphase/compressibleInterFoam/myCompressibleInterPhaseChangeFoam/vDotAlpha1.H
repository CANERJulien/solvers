Info<< "Calculating vDotAlphal\n" << endl;
dimensionedScalar zeroPressure
(
    "zeroPressure",
    dimPressure,
    0
);

volScalarField alpha1Coeff(1.0/rho1 - alpha1*(1.0/rho1 - 1.0/rho2));

/*
Pair<tmp<volScalarField>> mDotAlpha1
                            (
                                mcCoeff*max((p - pSat), zeroPress),
                                mvCoeff*min((p - pSat), zeroPress)
                            );

Pair<tmp<volScalarField>> vDotAlpha1 (alphalCoeff*mDotAlphal[0], alphalCoeff*mDotAlphal[1]);
*/
//volScalarField vDotcAlphal(alpha1Coeff*mcCoeff*max((p - pSat), zeroPressure));
volScalarField vDotcAlpha1
(
    alpha1Coeff*mcCoeff*sqrt(1/T)
   *mag(fvc::grad(alpha1))*unitLength
   *max((p - pSat), zeroPressure)
);
//volScalarField vDotvAlphal(alpha1Coeff*mvCoeff*min((p - pSat), zeroPressure));
volScalarField vDotvAlpha1
(
    alpha1Coeff*mvCoeff*sqrt(1/T)
   *mag(fvc::grad(alpha1))*unitLength
   *min((p - pSat), zeroPressure)
);

volScalarField vDotvmcAlpha1(vDotvAlpha1 - vDotcAlpha1);




Info<< "Volumetric condensation and vaporization rates \n"
    << "Min(vDotcAlpha1) = " << min(vDotcAlpha1).value()
    << " Max(vDotcAlpha1) = " << max(vDotcAlpha1).value() << "\n"
    << "Min(vDotvAlpha1) = " << min(vDotvAlpha1).value()
    << " Max(vDotvAlpha1) = " << max(vDotvAlpha1).value() << "\n"
    << "Min(vDotvmcAlpha1) = " << min(vDotvmcAlpha1).value()
    << " Max(vDotvmcAlpha1) = " << max(vDotvmcAlpha1).value()
    << endl;
/*
if (mesh.time().outputTime())
{
    volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));
    volScalarField vDotv
    (
    "vDotv",
    vDotvAlpha1*limitedAlpha1
    );
    volScalarField vDotc
    (
    "vDotc",
    vDotcAlpha1*(1-limitedAlpha1)
    );



    vDotv.write();
    vDotc.write();
}
*/
