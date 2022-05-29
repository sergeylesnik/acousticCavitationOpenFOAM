// * * * * * * * * * * * * * * * From 0/ Directory * * * * * * * * * * * * * //

    Info<< "Reading field Real Part of Acoustic Pressure\n" << endl;
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

    Info<< "Reading field Imaginary Part of Acoustic Pressure\n" << endl;
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

    Info<< "Defining field Absolute Acoustic Pressure\n" << endl;
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
        sqrt(sqr(PAcRe) + sqr(PAcIm))
    );

    Info<< "Reading field Phase of Acoustic Pressure\n" << endl;
    volScalarField phiPAc
    (
        IOobject
        (
            "phiPAc",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Acoustic Pressure Gradient\n" << endl;
    volVectorField G
    (
        IOobject
        (
            "G",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        G
    );
    psi.dimensions().reset(dimensionSet(0, 0, 0, 0, 0));

    // * * * * * * * * * * * * * * * Flow field  * * * * * * * * * * * * * //


    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
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

#   include "createPhi.H"

    scalar pRefValue = readScalar
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("pRefValue")
    );
    label pRefCell = readScalar
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("pRefCell")
    );

    mesh.schemesDict().setFluxRequired(p.name());

    // Transport model and viscosity model are created and parameters are
    // read from transportProperties dict.
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    Info<< "Reading field rAU if present\n" << endl;
    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    );

    Info << "Reading g\n" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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

    Info<< "Reading a constant Density\n" << endl;
    dimensionedScalar rho
    (
        transportProperties.lookup("rho")
    );

    Info<< "Creating a field Density since kinematicCloud assumes"
        << " compressible flow\n" << endl;
    volScalarField rhoField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rho
    );

    Info<< "Creating a field dynamic viscosity since kinematicCloud assumes"
        << " non-laminar flow\n" << endl;
    volScalarField muField("mu", turbulence->nu()*rho);


    // * * * * * * * * * * * * Fields for Damping calc * * * * * * * * * * //

    Info<< "Reading field Imaginary Wave Number\n" << endl;
    volScalarField kSqrIm
    (
        IOobject
        (
            "kSqrIm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field kSqrRe: real part of the squared wave number\n"
        << endl;
    volScalarField kSqrRe
    (
        IOobject
        (
            "kSqrRe",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        kSqrIm
    );

    Info<< "Creating field voidFrac (bubble void fraction)\n" << endl;
    volScalarField voidFrac
    (
        IOobject
        (
            "voidFrac",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField dKidPr
    (
        IOobject
        (
            "dKidPr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/(dimPressure*dimLength*dimLength)
    );

    // Construct others as copies resetting names
    volScalarField dKrdPr("dKrdPr", dKidPr);
    volScalarField dKidPi("dKidPi", dKidPr);
    volScalarField dKrdPi("dKrdPi", dKidPr);

    volScalarField dPAcRe
    (
        IOobject
        (
            "dPAcRe",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        PAcRe
    );

    volScalarField dPAcIm("dPAcIm", dPAcRe);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    // Initial injection is in the constructor.
    AcCavitationCloud<acCavitationParcel> bubbleCloud
    (
        "bubbleCloud",
        rhoField,
        U,
        muField,
        PAc,
        phiPAc,
        G,
        psi,
        g
    );

    const dictionary& ctrlDict
    (
        mesh.solutionDict().subDict("newtonRaphson")
    );

    IOdictionary derivedProperties
    (
        IOobject
        (
            "derivedProperties",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dictionary("derivedProperties")
    );

    linearDampingWaveNumber kSqrLinDamp
    (
        transportProperties,
        bubbleCloud.constProps().dict(),
        derivedProperties
    );

    derivedProperties.regIOobject::write();
