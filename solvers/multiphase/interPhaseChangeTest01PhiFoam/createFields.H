Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

//Add "T" dictionary entry here...
Info<< "Reading field T\n" << endl;

volScalarField T
(    
    IOobject
    (
       "T",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
     ),
     mesh
);

surfaceScalarField phirAlpha
            (
                IOobject
                (
                    "phirAlpha",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                    
                ),
               phi
                
);

//Done...



Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);


scalar xmin=-1e30;
scalar ymin=1e30;
scalar zmax=-1e30;

/*forAll (alpha1,cellI)
       {
            if ((U.mesh().C()[cellI].x()<=xmin) || (U.mesh().C()[cellI].y()<=ymin) || (U.mesh().C()[cellI].z()>=zmax))
		      {
			   xmin=U.mesh().C()[cellI].x();
			   ymin=U.mesh().C()[cellI].y();
			   zmax=U.mesh().C()[cellI].z();
			   pRefCell=cellI;
		      }
	    }*/
forAll (alpha1,cellI)
       {
            if (((U.mesh().C()[cellI].y()<=ymin) || (U.mesh().C()[cellI].z()>=zmax)) && ((U.mesh().C()[cellI].x()<=scalar(0))))
		      {
			   xmin=U.mesh().C()[cellI].x();
			   ymin=U.mesh().C()[cellI].y();
			   zmax=U.mesh().C()[cellI].z();
               //xmin=scalar(0);
               pRefCell=cellI; 
			  }
		    
	    }

p[pRefCell]=pRefValue;
Info<< "pRefCell\n" << U.mesh().C()[pRefCell].x() << "**" << U.mesh().C()[pRefCell].y() <<"**"  << U.mesh().C()[pRefCell].z() << "**"<< pRefValue << endl;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"

//Add "thermoPhyProperties" dictionary entry here...
Info<< "Reading thermoPhysProperties\n" << endl;

IOdictionary thermoPhysProperties
(    
    IOobject
    (
       "thermoPhysProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ,
       IOobject::NO_WRITE
     )
);


dimensionedScalar k1
(
     thermoPhysProperties.lookup("k1")
);

dimensionedScalar k2
(
     thermoPhysProperties.lookup("k2")
);

dimensionedScalar rho1cp1
(
     thermoPhysProperties.lookup("rho1cp1")
);

dimensionedScalar rho2cp2
(
     thermoPhysProperties.lookup("rho2cp2")
);

dimensionedScalar P0
(
     thermoPhysProperties.lookup("P0")
);

dimensionedScalar T0
(
     thermoPhysProperties.lookup("T0")
);

dimensionedScalar M
(
     thermoPhysProperties.lookup("M")
);

dimensionedScalar lv
(
     thermoPhysProperties.lookup("lv")
);

dimensionedScalar fc
(
     thermoPhysProperties.lookup("fc")
);

dimensionedScalar fv
(
     thermoPhysProperties.lookup("fv")
);


dimensionedScalar unitLength
(
     thermoPhysProperties.lookup("unitLength")
);

dimensionedScalar Kc
(
     thermoPhysProperties.lookup("Kc")
);

dimensionedScalar Rgas
(
"Rgas",
dimensionSet(1, 2, -2, -1, -1, 0, 0),
scalar(8.31447)
);



volScalarField rhocp
(    
    IOobject
    (
       "rhocp",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
     rho1cp1*alpha1+rho2cp2*alpha2
);

volScalarField k
(    
    IOobject
    (
       "k",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
    k1*min(alpha1,scalar(1))+k2*max(alpha2,scalar(0))
);


 volScalarField pSatu
(
    IOobject
    (
        "pSatu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
    mesh,
 P0
);

dimensionedScalar Q0
(
"Q0",
dimensionSet(1,-1,-3,0,0,0,0),
scalar(0.0)
); 

volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     mesh,
    Q0
 //(mDotAlphal[0]+mDotAlphal[1])*lv*mag(fvc::grad(alpha1))*unitLength
);




/*
pSatu = P0*exp((M*lv/Rgas)*(scalar(1)/T0-scalar(1)/T));
*/

/*
volScalarField mc
(    
    IOobject
    (
       "mc",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
     ),
     
     
    scalar(2)*fc/(scalar(2)-fc)*sqrt(M/(scalar(2*3.141459265358932384*8.31447)*T))*(P0*M*lv)/(scalar(8.31447)*T*T)*max(Tsatu-T,scalar(0))
);

volScalarField mv
(    
    IOobject
    (
       "mv",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
     ),
    scalar(2)*fv/(scalar(2)-fv)*sqrt(M/(scalar(2*3.141459265358932384*8.31447)*T))*(P0*M*lv)/(scalar(8.31447)*T*T)*min(Tsatu-T,scalar(0))
);

volScalarField m
(
    IOobject
    (
        "m",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mc+mv
);
*/








//Done...

