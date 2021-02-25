/*{
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
     ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
}*/

{
    e1 = mixture.thermo1().he();
    e2 = mixture.thermo2().he();
    
    e = (alpha1*rho1*e1 + alpha2*rho2*e2)/(alpha1*rho1 + alpha2*rho2);
 
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvm::div(rhoPhi, e)
      - fvm::laplacian(turbulence.alphaEff(), e)
      + divUp
      + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
      - (vDotvP - vDotcP)*(p - pSat)*deltaHvap // phase change source term
    );
 

    /*
    fvScalarMatrix EEqn1
    (
        fvm::ddt(alpha1*rho1, e1) + fvm::div(alpha1*rho1*U, e1)
      - fvm::laplacian(turbulence.alphaEff(), e1) // div(q)
      // phase change source term
    );

    fvScalarMatrix EEqn2
    (
        fvm::ddt(alpha2*rho2, e2) + fvm::div(alpha2*rho2*U, e2)
      - fvm::laplacian(turbulence.alphaEff(), e2) // div(q)
      // phase change source term
    );

   
   fvScalarMatrix EEqnMec
   (
        fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
      + divUp
   );
   */
    
     
    EEqn.relax();

    EEqn.solve();


}
