Info<< "Reading field PAcRe: real part of the acoustic pressure amplitude\n"
    << endl;
volScalarField PAcRe
(
    IOobject
    (
        "PAcRe",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field PAcIm: imaginary part of the acoustic pressure"
    << "amplitude\n" << endl;
volScalarField PAcIm
(
    IOobject
    (
        "PAcIm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field PAc: absolute acoustic pressure amplitude\n" << endl;
volScalarField PAc
(
    IOobject
    (
        "PAc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(pow(PAcRe, 2) + pow(PAcIm, 2))
);

Info<< "Creating working coupled solution field blockT\n" << endl;
volVector2Field blockT
(
    IOobject
    (
        "blockT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector2("zero", dimless, vector2::zero)
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading a constant kSqrRe: real part of squared wave"
    << "number according to the linear damping theory\n" << endl;
dimensionedScalar kSqrRe
(
    transportProperties.lookup("kSqrRe")
);

Info<< "Reading a constant kSqrIm: imaginary part of squared wave"
    << "number according to the linear damping theory\n" << endl;
dimensionedScalar kSqrIm
(
    transportProperties.lookup("kSqrIm")
);

Info<< "Reading MUMPSSettings\n" << endl;
IOdictionary MUMPSSettings
(
    IOobject
    (
     "MUMPSSettings",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ
    )
);
