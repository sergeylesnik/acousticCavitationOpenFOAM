// The constructor has to be included in header otherwise linker error.
inline Foam::acCavitationParcel::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    Foam::KinematicParcel<acCavitationParcel>::constantProperties(parentDict),
    freq_(this->dict().lookup("freq")),
    pInf_(this->dict().lookup("pInf")),
    sigma_(this->dict().lookup("sigma")),
    DThermal_(this->dict().lookup("DThermal")),
    gamma_(this->dict().lookup("gamma")),
    spread_(readScalar(this->dict().lookup("cavitationInceptionSmoothing"))),
    nMaxParticlesPerParcel_
    (
        readScalar(this->dict().lookup("nMaxParticlesPerParcel"))
    ),
    nMinParcelsPerCell_
    (
        this->dict().lookupOrDefault<scalar>("nMinParcelsPerCell", 1.0)
    )
    {}


// The constructor has to be included in header otherwise linker error.
inline Foam::acCavitationParcel::trackData::trackData
(
    AcCavitationCloud<acCavitationParcel>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const interpolation<scalar>& PAcInterp,
    const interpolation<scalar>& argPAcInterp,
    const interpolation<vector>& GInterp,
    const interpolation<vector>& argGInterp,
    const vector& g
)
:
    KinematicParcel<acCavitationParcel>::trackData
    (
        cloud,
        constProps,
        rhoInterp,
        UInterp,
        muInterp,
        g
    ),
    cloud_(cloud),
    constProps_(constProps),
    ReAbsTol_
    (
        constProps.dict().lookupOrDefault<scalar>
        (
            "ReAbsTol", 10.0
        )
    ),
    PAcInterp_(PAcInterp),
    argPAcInterp_(argPAcInterp),
    GInterp_(GInterp),
    argGInterp_(argGInterp),
    RAvTable_
    (
        interpolation2DTable<scalar>
        (
            cloud.particleProperties().subDict("RAvTable")
        )
    ),
    VAvTable_
    (
        interpolation2DTable<scalar>
        (
            cloud.particleProperties().subDict("VAvTable")
        )
    ),
    IcTable_
    (
        interpolation2DTable<scalar>
        (
            cloud.particleProperties().subDict("IcTable")
        )
    ),
    IsTable_
    (
        interpolation2DTable<scalar>
        (
            cloud.particleProperties().subDict("IsTable")
        )
    ),
    patchRndInjectorPtr_(nullptr),
    parcelsAddedConti_(0),
    massAddedConti_(0.0),
    nShiftedParcels_(0),
    nParcelShifts_(0)
{
    try 
    {
        patchRndInjectorPtr_ = 
            dynamic_cast
            <PatchThresholdRndInjection<KinematicCloud<acCavitationParcel> >* >
            (& cloud.injection());
    }
    catch(std::bad_cast& e)
    {
        // do nothing since already a nullptr
    }
}


// Copy constructor
inline Foam::acCavitationParcel::trackData::trackData
(
    Foam::acCavitationParcel::trackData& td
)
:
    KinematicParcel<acCavitationParcel>::trackData
    (
        td.cloud(),
        td.constProps(),
        td.rhoInterp(),
        td.UInterp(),
        td.muInterp(),
        td.g()
    ),
    cloud_(td.cloud()),
    constProps_(td.constProps()),
    ReAbsTol_
    (
        constProps_.dict().lookupOrDefault<scalar>
        (
            "ReAbsTol", 10
        )
    ),
    PAcInterp_(td.PAcInterp()),
    argPAcInterp_(td.argPAcInterp()),
    GInterp_(td.GInterp()),
    argGInterp_(td.argGInterp()),
    RAvTable_(td.RAvTable()),
    VAvTable_(td.VAvTable()),
    IcTable_(td.IcTable()),
    IsTable_(td.IsTable()),
    patchRndInjectorPtr_(nullptr),
    parcelsAddedConti_(0),
    massAddedConti_(0.0),
    nShiftedParcels_(0),
    nParcelShifts_(0)
{
    try 
    {
        patchRndInjectorPtr_ = 
            dynamic_cast
            <PatchThresholdRndInjection<KinematicCloud<acCavitationParcel> >* >
            (& cloud_.injection());
    }
    catch(std::bad_cast& e)
    {
        // do nothing since already a nullptr
    }
}


// Used in the PatchThresholdRndInjection
inline Foam::acCavitationParcel::acCavitationParcel
(
    KinematicCloud<acCavitationParcel>& owner,
    const vector& position,
    const label cellI
)
:
    KinematicParcel<acCavitationParcel>(owner, position, cellI),
    dAv_(this->d_),
    VAv_(0.0),
    lastCell_(cellI),
    nSmallDeltaStepFractions_(0),
    FBjPri_(vector::zero),
    PAc_(0.0)
{}


inline Foam::acCavitationParcel::acCavitationParcel
(
    KinematicCloud<acCavitationParcel>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const constantProperties& constProps
)
:
    KinematicParcel<acCavitationParcel>
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        constProps
    ),
    dAv_(this->d_),
    VAv_(0.0),
    lastCell_(cellI),
    nSmallDeltaStepFractions_(0),
    FBjPri_(vector::zero),
    PAc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

inline const Foam::dimensionedScalar&
Foam::acCavitationParcel::constantProperties::sigma() const
{
    return sigma_;
}


inline const Foam::dimensionedScalar&
Foam::acCavitationParcel::constantProperties::pInf() const
{
    return pInf_;
}


inline Foam::scalar
Foam::acCavitationParcel::constantProperties::nMaxParticlesPerParcel() const
{
    return nMaxParticlesPerParcel_;
}


inline Foam::scalar
Foam::acCavitationParcel::constantProperties::nMinParcelsPerCell() const
{
    return nMinParcelsPerCell_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

inline Foam::AcCavitationCloud<Foam::acCavitationParcel>&
Foam::acCavitationParcel::trackData::cloud()
{
    return cloud_;
}


inline const Foam::AcCavitationCloud<Foam::acCavitationParcel>&
Foam::acCavitationParcel::trackData::cloud() const
{
    return cloud_;
}


inline const typename Foam::acCavitationParcel::constantProperties&
Foam::acCavitationParcel::trackData::constProps() const
{
    return constProps_;
}


inline Foam::scalar
Foam::acCavitationParcel::trackData::ReAbsTol() const
{
    return ReAbsTol_;
}


inline const Foam::interpolation<Foam::scalar>&
Foam::acCavitationParcel::trackData::PAcInterp() const
{
    return PAcInterp_;
}


inline const Foam::interpolation<Foam::scalar>&
Foam::acCavitationParcel::trackData::argPAcInterp() const
{
    return argPAcInterp_;
}


inline const Foam::interpolation<Foam::vector>&
Foam::acCavitationParcel::trackData::GInterp() const
{
    return GInterp_;
}


inline const Foam::interpolation<Foam::vector>&
Foam::acCavitationParcel::trackData::argGInterp() const
{
    return argGInterp_;
}


inline const Foam::interpolation2DTable<Foam::scalar>&
Foam::acCavitationParcel::trackData::RAvTable() const
{
    return RAvTable_;
}


inline const Foam::interpolation2DTable<Foam::scalar>&
Foam::acCavitationParcel::trackData::VAvTable() const
{
    return VAvTable_;
}


inline const Foam::interpolation2DTable<Foam::scalar>&
Foam::acCavitationParcel::trackData::IcTable() const
{
    return IcTable_;
}


inline const Foam::interpolation2DTable<Foam::scalar>&
Foam::acCavitationParcel::trackData::IsTable() const
{
    return IsTable_;
}


inline Foam::scalar
Foam::acCavitationParcel::trackData::VAvTableInterpolate
(
    scalar REqu,
    scalar PAc
) const
{
    return VAvTable_(REqu, PAc);
}


inline Foam::scalar
Foam::acCavitationParcel::trackData::IcTableInterpolate
(
    scalar REqu,
    scalar PAc
) const
{
    return IcTable_(REqu, PAc);
}


inline Foam::scalar
Foam::acCavitationParcel::trackData::IsTableInterpolate
(
    scalar REqu,
    scalar PAc
) const
{
    return IsTable_(REqu, PAc);
}


inline Foam::PatchThresholdRndInjection
    <
        Foam::KinematicCloud<Foam::acCavitationParcel> 
    >*
Foam::acCavitationParcel::trackData::patchRndInjectorPtr()
{
    return patchRndInjectorPtr_;
}


inline Foam::label&
Foam::acCavitationParcel::trackData::parcelsAddedConti()
{
    return parcelsAddedConti_;
}


inline Foam::scalar&
Foam::acCavitationParcel::trackData::massAddedConti()
{
    return massAddedConti_;
}


inline Foam::label&
Foam::acCavitationParcel::trackData::nShiftedParcels()
{
    return nShiftedParcels_;
}


inline Foam::label&
Foam::acCavitationParcel::trackData::nParcelShifts()
{
    return nParcelShifts_;
}


// * * * * * * * * acCavitationParcel Member Functions * * * * * * * * * * * //

inline Foam::scalar
Foam::acCavitationParcel::REqu() const
{
    return d_/2;
}


inline Foam::scalar
Foam::acCavitationParcel::dAv() const
{
    return dAv_;
}


inline Foam::vector
Foam::acCavitationParcel::FBjPri() const
{
    return FBjPri_;
}


inline Foam::scalar
Foam::acCavitationParcel::PAc() const
{
    return PAc_;
}


inline Foam::label
Foam::acCavitationParcel::lastCell() const
{
    return lastCell_;
}


inline Foam::label&
Foam::acCavitationParcel::lastCell()
{
    return lastCell_;
}


inline Foam::label
Foam::acCavitationParcel::nSmallDeltaStepFractions() const
{
    return nSmallDeltaStepFractions_;
}