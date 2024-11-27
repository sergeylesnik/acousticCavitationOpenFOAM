#include "AcCavitationCloud.H"
#include "interpolationCellPoint.H"
#include "acCavitationParcel.H"
#include "pdf.H"
#include "meshTools.H"

#include "vectorList.H"
#include "PatchThresholdRndInjection.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::injectParcelInCell
(
    const fvMesh& mesh,
    const vector& deltaCoord,
    const vector& minCoord,
    const meshSearch& searchEngine,
    const autoPtr<pdf>& parcelPdf,
    const label cellI,
    const scalar VCellI,
    scalar nP,
    label& injCounter,
    scalar& voidFracCellI
)
{
    point injPosition;

    // Since the cuboid is larger than the cell, the point might not lie
    // inside the cell. In this case randomize the position again.
    do
    {
        // Inject position is the sum of the lower conner of the cuboid
        // and a random position inside the cuboid.
        injPosition =
            minCoord
            + cmptMultiply(this->rndGen().vector01(), deltaCoord);
    }
    while (!searchEngine.pointInCell(injPosition, cellI));

    // In case of wedge, constrain the circumferential direction of the
    // injection position.
    meshTools::constrainDirection(mesh, mesh.solutionD(), injPosition);

    ParcelType* pPtr = new ParcelType(*this, injPosition, cellI);
    pPtr->d() = parcelPdf->sample();
    pPtr->nParticle() = (nP < 1.0) ? 1.0 : nP;

    // Set other parcel properties to default
    checkParcelProperties(*pPtr, 0.0, false);
    this->addParticle(pPtr);

    injCounter++;
    voidFracCellI += pPtr->volume() * pPtr->nParticle() / VCellI;

}

template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::preEvolve()
{
    this->dispersion().cacheFields(true);
    acCavForces_.cacheFields(true);
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::evolveCloud()
{
    autoPtr<interpolation<scalar> > rhoInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > PAcInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        PAc()
    );

    autoPtr<interpolation<scalar> > argPAcInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        argPAc()
    );

    autoPtr<interpolation<vector> > GInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        G()
    );

    autoPtr<interpolation<vector> > argGInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        argG()
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp(),
        UInterp(),
        muInterp(),
        PAcInterp(),
        argPAcInterp(),
        GInterp(),
        argGInterp(),
        this->g().value()
    );

    this->injection().inject(td);

    if (this->coupled())
    {
        resetSourceTerms();
    }

    Cloud<ParcelType>::move(td);

    info(td);

    parcelsAddedConti_ = td.parcelsAddedConti();
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::postEvolve()
{
    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);
    acCavForces_.cacheFields(false);

    this->postProcessing().post();

    // Update number of parcels per cell
    parcelsPerCell_.field() = 0.0;
    labelList nppc = nParcelsPerCell();
    forAll(parcelsPerCell_, cellI)
    {
       parcelsPerCell_[cellI] = nppc[cellI];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::AcCavitationCloud<ParcelType>::AcCavitationCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const volScalarField& PAc,
    const volScalarField& argPAc,
    const volVectorField& G,
    const volVectorField& argG,
    const dimensionedVector& g,
    bool readFields
)
:
    KinematicCloud<ParcelType>
    (
        cloudName,
        rho,
        U,
        mu,
        g,
        false
    ),
    constProps_(this->particleProperties()),
    PAc_(PAc),
    argPAc_(argPAc),
    G_(G),
    argG_(argG),
    dampingViscous_
    (
        this->particleProperties().subDict("damping").lookup("viscous")
    ),
    dampingThermal_
    (
        this->particleProperties().subDict("damping").lookup("thermal")
    ),
    dampingRadiation_
    (
        this->particleProperties().subDict("damping").lookup("radiation")
    ),
    injectDomainInitDict_
    (
        this->particleProperties().subDict("domainInitialInjection")
    ),
    injectDomainInit_
    (
        injectDomainInitDict_.lookup("active")
    ),
    injectDomainFlagDict_
    (
        IOobject
        (
            "domainInitInjectFlag",
            this->db().time().timeName(),
            "uniform"/cloud::prefix/this->name(),
            this->db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    PiViscous_
    (
        IOobject
        (
            this->name() + "PiViscous",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    PiThermal_
    (
        IOobject
        (
            this->name() + "PiThermal",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    PiRadiation_
    (
        IOobject
        (
            this->name() + "PiRadiation",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    PiTotal_
    (
        IOobject
        (
            this->name() + "PiTotal",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    parcelsPerCell_
    (
        IOobject
        (
            this->name() + "ParcelsPerCell",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    PiViTablePtr_(nullptr),
    PiThTablePtr_(nullptr),
    acCavForces_
    (
        this->mesh(),
        this->particleProperties(),
        g.value()
    ),
    parcelsAddedConti_(0)
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }

    if (dampingViscous_)
    {
        PiViTablePtr_.set
        (
            new interpolation2DTable<scalar>
            (
                this->particleProperties().subDict("PiViscousTable")
            )
        );
    }

    if (dampingThermal_)
    {
        PiThTablePtr_.set
        (
            new interpolation2DTable<scalar>
            (
                this->particleProperties().subDict("PiThermalTable")
            )
        );
    }

    injectDomainInit();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::AcCavitationCloud<ParcelType>::~AcCavitationCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    KinematicCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    if (!fullyDescribed)
    {
        parcel.lastCell() = parcel.cell();
    }
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::resetSourceTerms()
{
    KinematicCloud<ParcelType>::resetSourceTerms();
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::evolve()
{
    if (this->active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();
    }
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::injectDomainInit()
{
    Info<< "Performing initialization of the bubbles in the domain" << endl;

    // Check the inject switches.
    word flagName("domainInjectionAlreadyDone");
    Switch injectDomainFlag
    (
        injectDomainFlagDict_.lookupOrDefault
        (
            flagName, false
        )
    );

    // Abort injection if any of the switches is set accordingly.
    if (!injectDomainInit_)
    {
        return;
    }
    else if (injectDomainFlag)
    {
        WarningIn("AcCavitationCloud<ParcelType>::injectDomainInit()")
            << "Flag identifying that the injection at the initialization"
            << " has already been done is set and no init injection will"
            << " be performed. If this is not the desired behaviour, change"
            << " the flag under " << nl << "    "
            << injectDomainFlagDict_.instance()
               /injectDomainFlagDict_.local()
               /injectDomainFlagDict_.name()
            << endl;

        return;
    }

    dictionary dictToRead;
    if (Switch(injectDomainInitDict_.lookup("useSettingsFromInjectionModel")))
    {
        dictToRead = this->injection().coeffDict();
        if (dictToRead == dictionary::null)
        {
            FatalErrorIn("AcCavitationCloud<ParcelType>::injectDomainInit()")
                << "useSettingsFromInjectionModel is on in "
                << injectDomainInitDict_.name()
                << " but injection model has no coeff dictionary"
                << abort(FatalError);
        }
    }
    else
    {
        dictToRead = injectDomainInitDict_;
    }

    const scalar voidFracInit =
        readScalar(dictToRead.lookup("thresholdVoidFrac"));

    const autoPtr<pdf> parcelInitPdf =
        pdf::New(dictToRead.subDict("parcelPDF"), this->rndGen());

    tmp<volScalarField> voidFrac = this->theta();
    const fvMesh& mesh = this->mesh();
    meshSearch searchEngine(mesh);
    label injCounter = 0;

    // Minimal possible volume of a bubble
    scalar minParticleVolume =
        mathematicalConstant::pi / 6.0 * pow3(parcelInitPdf->minValue());

    // Find out PDFs median diameter and corresponding volume by sampling the
    // PDF
    scalar dMedian = 0.0;
    label sampleSize = 1000;
    for (int i=0; i < sampleSize; i++)
    {
        dMedian += parcelInitPdf->sample();
    }
    dMedian /= sampleSize;
    Info<< "Median diameter of bubble PDF is " << dMedian << endl;
    scalar medianParticleVolume = mathematicalConstant::pi/6.0*pow3(dMedian);

    scalar nMaxParticlesPerParcel = this->constProps().nMaxParticlesPerParcel();
    scalar nMinParcelsPerCell = this->constProps().nMinParcelsPerCell();

    bool isMonodisperse =
        (parcelInitPdf->minValue() == parcelInitPdf->maxValue());
    if (isMonodisperse)
    {
        Info<< "AcCavitationCloud: monodisperse case; nParticle of last"
            << " parcel in a cell will be adjusted to fit threshold on average"
            << endl;
    }

    forAll(mesh.C(), cellI)
    {
        const scalar VCellI = mesh.cellVolumes()[cellI];
        scalar voidFracMin = minParticleVolume/VCellI;
        scalar voidFracMedian = medianParticleVolume/VCellI;

        // Proceed with the next cell if the minimal possible void fraction is
        // greater than the given threshold
        if (voidFracMin > voidFracInit)
        {
            continue;
        }

        // Find out the min and max coordinates of all the cell vertices.
        // By doing this we virtually put the cell into a hexahedron
        // (rectangular cuboid)
        const labelList& cellPointsInd = mesh.cellPoints()[cellI];
        List<point> cellPoints;
        forAll(cellPointsInd, i)
        {
            point pt = mesh.points()[cellPointsInd[i]];
            cellPoints.append(pt);
        }
        vector minCoord = min(cellPoints);
        vector deltaCoord = max(cellPoints);
        deltaCoord -= minCoord;

        scalar& voidFracCellI = voidFrac()[cellI];

        // Calculate such nParticle for each parcel that at least one parcel of // the median size fits into each cell.
        label nP = round(voidFracInit / (voidFracMedian*nMinParcelsPerCell));
        if (nP < 1.0)
        {
            nP = 1.0;
        }
        else if (nP > nMaxParticlesPerParcel)
        {
            nP = nMaxParticlesPerParcel;
        }

        // Fill up the cell with new parcels until the threshold void fraction
        // is reached. Assume that the threshold/cell volumes are large enough
        // and a bubble of random size from the pdf may be injected.
        while (voidFracCellI < voidFracInit - nP*voidFracMedian)
        {

            // Randomize number of particles in the upper range of nP in order
            // to avoid cell checkerboard pattern of void fraction.
            label nPCorr = this->rndGen().integer(round(0.9*nP), nP);

            injectParcelInCell
            (
                mesh,
                deltaCoord,
                minCoord,
                searchEngine,
                parcelInitPdf,
                cellI,
                VCellI,
                nPCorr,
                injCounter,
                voidFracCellI
            );
        }

        if (isMonodisperse)
        {
            // Compute nP to fit the threshold void fraction properly.
            scalar pVolume =
                mathematicalConstant::pi/6.0*pow3(parcelInitPdf->minValue());
            nP = round((voidFracInit - voidFracCellI)*VCellI/pVolume);

            if (nP > 1.0)
            {
                label nPCorr = nP;
                injectParcelInCell
                (
                    mesh,
                    deltaCoord,
                    minCoord,
                    searchEngine,
                    parcelInitPdf,
                    cellI,
                    VCellI,
                    nPCorr,
                    injCounter,
                    voidFracCellI
                );
            }
        }
    }

    voidFrac.clear();

    // Set the init injection flag to true.
    injectDomainFlagDict_.set(flagName, true);

    injCounter = returnReduce(injCounter, sumOp<label>());
    Info<< "Total of " << injCounter << " parcels have been injected" << endl;

}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    updateMesh();
}

template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::updateMesh()
{
    KinematicCloud<ParcelType>::updateMesh();

    forAllIter(typename Cloud<ParcelType>, *this, pIter)
    {
        pIter().lastCell() = pIter().cell();
    }
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::info
(
    typename ParcelType::trackData& td
) const
{
    Info<< "Cloud: " << this->name() << nl
        << "    Number of parcels added by injector         = "
        << this->injection().parcelsAddedTotal() << nl
        << "    Mass introduced by injector                 = "
        << this->injection().massInjected() << nl
        << "    Number of parcels added continuously        = "
        << returnReduce(td.parcelsAddedConti(), sumOp<label>()) << nl
        << "    Mass introduced continuously                = "
        << returnReduce(td.massAddedConti(), sumOp<scalar>()) << nl
        << "    Current number of parcels                   = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system                      = "
        << returnReduce(this->massInSystem(), sumOp<scalar>()) << nl
        << "    Number of parcels shifted at least once     = "
        << returnReduce(td.nShiftedParcels(), sumOp<label>()) << nl
        << "    Number of parcel shifts towards cell center = "
        << returnReduce(td.nParcelShifts(), sumOp<label>()) << nl
        << endl;
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::computeImagWaveNumber
(
    DimensionedField<scalar, volMesh>& kSqrIm,
    DimensionedField<scalar, volMesh>& kSqrRe,
    const DimensionedField<scalar, volMesh>& PAc,
    const DimensionedField<scalar, volMesh>& voidFrac,
    linearDampingWaveNumber& ldwn,
    stepSmoother& ss
)
{
    resetAcousticFields();

    // Calculation of the real part of the squared wave number cannot be split
    // over individual bubbles like the imaginary one. Thus, compute an
    // average equilibrium radius and compute Blake threshold from this.
    scalar REquAv = 0.0;
    label nPTot = 0;

    forAllConstIter(typename AcCavitationCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();
        const scalar REqu = p.REqu();
        const scalar nP = p.nParticle();
        const scalar PAcCell = PAc[cellI];

        scalar PBlake = ldwn.computeBlakePressure(REqu);
        scalar nuclCoeff = ss.smoothStep(PBlake, PAcCell);

        if (nuclCoeff > SMALL)
        {
            if (dampingViscous_)
            {
                PiViscous_[cellI] +=
                    nP*nuclCoeff*PiViTablePtr_()(REqu, PAcCell);
            }
            if (dampingThermal_)
            {
                PiThermal_[cellI] +=
                    nP*nuclCoeff*PiThTablePtr_()(REqu, PAcCell);
            }
        }

        REquAv += REqu*nP;
        nPTot++;
    }

    if (dampingViscous_)
    {
        PiTotal_ += PiViscous_;
    }
    if (dampingThermal_)
    {
        PiTotal_ += PiThermal_;
    }

    // Make it the same on all procs
    REquAv = returnReduce(REquAv, sumOp<scalar>());
    scalar gSize = returnReduce(this->size(), sumOp<scalar>());
    REquAv /= (nPTot*gSize + SMALL);
    REquAv = (REquAv < SMALL) ? SMALL : REquAv;

    // Imaginary part
    const dimensionedScalar& omega = ldwn.omegaDim();
    const dimensionedScalar& rhoc = ldwn.rhoLiquidDim();
    const DimensionedField<scalar, volMesh>& cellVolumes = this->mesh().V();

    dimensionedScalar PAcUnit("PAcUnit", dimPressure, 1.0);
    kSqrIm =
        -2.0*rhoc*omega*PiTotal_.dimensionedInternalField()
      / (sqr(PAc + PAcUnit*SMALL)*cellVolumes);

    // Real part
    scalar PBlakeAv = ldwn.computeBlakePressure(REquAv);
    scalar kSqrReNoCav = ldwn.kSqrReNoCav();
    ldwn.computeDampingAndFreq(REquAv);
    forAll(kSqrRe, cellI)
    {
        scalar nuclCoeff = ss.smoothStep(PBlakeAv, PAc[cellI]);
        scalar kSqrReCav = ldwn.realPartWithVoidFrac(voidFrac[cellI]);
        kSqrRe[cellI] = (1.0 - nuclCoeff)*kSqrReNoCav + nuclCoeff*kSqrReCav;
    }

//label I = 1;
//Info<< "PiViscous_[I] = " << PiViscous_[I] << nl
//    << "testCellIndex = " << I << nl
//    << "PiViTable[I] = " << PiViTablePtr_()(2e-6, PAc[I]) << nl
//    << "parcelsPerCell_[I] = " << parcelsPerCell_[I] << nl
//    << "voidFrac[I] = "
//    << 3.1415/6.0*pow3(4e-6)*parcelsPerCell_[I]/this->mesh().V()[I] << nl
//    << "PAc[I] = " << PAc[I] << nl
//    << "rhoc = " << rhoc << nl
//    << "omega = " << omega << nl
//    << "this->mesh().V() = " << this->mesh().V()[I] << nl
//    << "kSqrIm[I] = " << kSqrIm[I] << nl << endl;
}


template<class ParcelType>
void Foam::AcCavitationCloud<ParcelType>::resetAcousticFields()
{
    PiViscous_.field() = 0.0;
    PiThermal_.field() = 0.0;
    PiRadiation_.field() = 0.0;
    PiTotal_.field() = 0.0;
}


template<class ParticleType>
Foam::labelList Foam::AcCavitationCloud<ParticleType>::nParcelsPerCell() const
{
    labelList nppc(this->pMesh().nCells(), 0);

    forAllConstIter(typename AcCavitationCloud<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        const label cellI = p.cell();

        // Check
        if (cellI < 0 || cellI >= this->pMesh().nCells())
        {
            FatalErrorIn
            (
                "Foam::AcCavitationCloud<ParticleType>::nParcelsPerCell()"
            )
                << "Illegal cell number " << cellI
                << " at position " << p.position() << nl
                << "Cell number should be between 0 and "
                << this->pMesh().nCells()-1 << nl
                << "On this mesh the particle should be in cell "
                << this->pMesh().findCell(p.position())
                << exit(FatalError);
        }

        nppc[cellI]++;
    }

    return nppc;
}

// ************************************************************************* //
