/* Essayer de modifier Su et Sp en y ajoutant vDotAlpha1 */

Info<< "Calculating vDotAlpha1\n" << endl;
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

Pair<tmp<volScalarField>> vDotAlpha1 (alpha1Coeff*mDotAlpha1[0], alpha1Coeff*mDotAlpha1[1]);
*/
//volScalarField vDotcAlpha1(alpha1Coeff*mcCoeff*max((p - pSat), zeroPressure));
volScalarField vDotcAlpha1(alpha1Coeff*mcCoeff*sqrt(1/T)*max((p - pSat), zeroPressure));
//volScalarField vDotvAlpha1(alpha1Coeff*mvCoeff*min((p - pSat), zeroPressure));
volScalarField vDotvAlpha1(alpha1Coeff*mvCoeff*sqrt(1/T)*min((p - pSat), zeroPressure));
volScalarField vDotvmcAlpha1(vDotvAlpha1 - vDotcAlpha1);

forAll(alpha1Coeff, celli)
{
    Sp[celli] += vDotvmcAlpha1[celli];
    Su[celli] += vDotcAlpha1[celli];
}



Info<< "Volumetric condensation and vaporization rates \n"
    << "Min(vDotcAlpha1) = " << min(vDotcAlpha1).value()
    << " Max(vDotcAlpha1) = " << max(vDotcAlpha1).value() << "\n"
    << "Min(vDotvAlpha1) = " << min(vDotvAlpha1).value()
    << " Max(vDotvAlpha1) = " << max(vDotvAlpha1).value() << "\n"
    << "Min(vDotvmcAlpha1) = " << min(vDotvmcAlpha1).value()
    << " Max(vDotvmcAlpha1) = " << max(vDotvmcAlpha1).value()
    << endl;


