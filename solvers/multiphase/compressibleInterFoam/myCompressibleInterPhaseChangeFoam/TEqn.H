#include "QT.H"
{
/*
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp()// - contErr/rho*p
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
          +  fvm::Sp(vDotvT, T) - vDotvT*TSat // phasechange source term
          +  fvm::Sp(vDotcT, T) - vDotcT*TSat
     ==
        fvOptions(rho, T)
    );
*/

    volScalarField rhoCp
    (
        alpha1*rho1*mixture.thermo1().Cp()
      + (scalar(1)-alpha1)*rho2*mixture.thermo2().Cp()  
    );

    surfaceScalarField rhoCpPhi
    (
        fvc::interpolate(rhoCp)*phi
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T) + fvm::div(rhoCpPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp()// - contErr/rho*p
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K
        )
      + fvm::Sp(QT, T) -QT*TSat //phase change source term
     ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
}
