
Info<< "Reading parameters for Clapeyron and Hertz-Knudsen\n" << endl;

IOdictionary phaseChangeParameters
(
    IOobject
    (
        "phaseChangeParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Clapeyron parameters
dimensionedScalar P0
(
    "P0",
    dimPressure,
    //dimensionSet(1,-1,-2,0,0,0,0),
    //phaseChangeParameters.subDict("Clapeyron").lookup("P0")
    phaseChangeParameters.subDict("Clapeyron")
);

dimensionedScalar T0
(
    "T0",
    dimTemperature,
    //phaseChangeParameters.subDict("Clapeyron").lookup("T0")
    phaseChangeParameters.subDict("Clapeyron")
);

dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    phaseChangeParameters.subDict("Clapeyron")
);

dimensionedScalar deltaHvap
(
    "deltaHvap",
    dimEnergy/dimMass,
    phaseChangeParameters.subDict("Clapeyron")
);

dimensionedScalar W
(
    "W",
    dimMass/dimMoles,
    phaseChangeParameters.subDict("Clapeyron")
);

// Hertz-Knudsen parameters
dimensionedScalar fc
(
    "fc",
    dimless,
    phaseChangeParameters.subDict("Hertz_Knudsen")
);
dimensionedScalar fv
(
    "fv",
    dimless,
    phaseChangeParameters.subDict("Hertz_Knudsen")
);
dimensionedScalar unitLength
(
    "unitLength",
    dimLength,
    phaseChangeParameters.subDict("Hertz_Knudsen")
);


dimensionedScalar mcCoeff
(
    //scalar(2)*fc/(scalar(2)-fc)*sqrt(W/(scalar(2*3.141459)*R*T0))/unitLength
    //scalar(2)*fc/(scalar(2)-fc)*sqrt(W/(scalar(2*3.141459)*R))/unitLength //use field T instead of T0
    scalar(2)*fc/(scalar(2)-fc)*sqrt(W/(constant::mathematical::twoPi*R))/unitLength //use field T instead of T0
);
dimensionedScalar mvCoeff
(
    //scalar(2)*fv/(scalar(2)-fv)*sqrt(W/(scalar(2*3.141459)*R*T0))/unitLength
    //scalar(2)*fv/(scalar(2)-fv)*sqrt(W/(scalar(2*3.141459)*R))/unitLength
    scalar(2)*fv/(scalar(2)-fv)*sqrt(W/(constant::mathematical::twoPi*R))/unitLength
);
