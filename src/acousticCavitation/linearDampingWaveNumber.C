#include "linearDampingWaveNumber.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::complex<Foam::scalar> Foam::linearDampingWaveNumber::i_(0.0, 1.0);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::linearDampingWaveNumber::computeComplex
(
    const Foam::scalar voidFrac,
    const Foam::scalar REqu
)
{            
    computeDampingAndFreq(REqu);

    // Complex wave number
    std::complex<scalar> kSqr =
        sqr(omega_/cSound_)
      + 3.0*voidFrac*sqr(omega_)
      / (sqr(REqu)*(omega0Sqr_ - sqr(omega_) + 2.0*i_*b_*omega_));

    kSqrReCav_ = kSqr.real();
    kSqrIm_ = kSqr.imag();

    // Speed of sound due to linear damping
    std::complex<scalar> cSoundLinDamp = omega_/std::sqrt(kSqr);

    // Put into Foam format to save it in dictionary
    Foam::complex kSqrFoam(kSqr.real(), kSqr.imag());
    Foam::complex cSoundLinDampFoam
    (
        cSoundLinDamp.real(),
        cSoundLinDamp.imag()
    );

    // Assume that this is done only in case of constant void fraction and it
    // is computed once. Thus, save value to the dictionary.
    derivedProps_.add("kSqrCmpxLinearDamp", kSqrFoam);
    derivedProps_.add("cSoundCmpxLinearDamp", cSoundLinDampFoam);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearDampingWaveNumber::linearDampingWaveNumber
(
    const dictionary& transportProps,
    const dictionary& bubbleProps,
    dictionary& derivedProps,
    dimensionedScalar REqu,
    dimensionedScalar voidFrac
)
:
    pi_(mathematicalConstant::pi),
    transportProps_(transportProps),
    bubbleProps_(bubbleProps),  
    derivedProps_(derivedProps), 
    freqDim_(bubbleProps_.lookup("freq")),
    pInfDim_(bubbleProps_.lookup("pInf")),
    REquDim_(REqu),
    voidFracDim_(voidFrac),
    rhocDim_(transportProps_.lookup("rho")),
    cSoundDim_(transportProps_.lookup("cSound")),
    nuDim_(transportProps_.lookup("nu")),
    sigmaDim_(bubbleProps_.lookup("sigma")),
    DThermalDim_(bubbleProps_.lookup("DThermal")),
    gammaDim_(bubbleProps_.lookup("gamma")),
    omegaDim_(2.0*pi_*freqDim_),
    muDim_(rhocDim_*nuDim_),
    freq_(freqDim_.value()),
    pInf_(pInfDim_.value()),
    REqu_(REquDim_.value()),
    voidFrac_(voidFracDim_.value()),
    rhoc_(rhocDim_.value()),
    cSound_(cSoundDim_.value()),
    nu_(nuDim_.value()),
    sigma_(sigmaDim_.value()),
    DThermal_(DThermalDim_.value()),
    gamma_(gammaDim_.value()),
    omega_(omegaDim_.value()),
    mu_(muDim_.value()),
    kSqrReNoCav_(sqr(omega_/cSound_))
{
    derivedProps_.add("omega", omegaDim_);
    
    if(voidFrac_ != 0.0 && REqu_ != 0.0)
    {
        compute();
    }
}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearDampingWaveNumber::compute()
{
    if (voidFrac_ == 0.0 || REqu_ == 0.0)
    {
        FatalErrorIn
        (
            "linearDampingWaveNumber::Compute()"
        )   << "Bubble void fraction was not set in the constructor call "
            << "or was set to 0. Wave number calculation would be "
            << "erroneous. You might want to use compute(scalar voidFrac)."
            << abort(FatalError);
    }

    computeDampingAndFreq(REqu_);
    computeComplex(voidFrac_, REqu_);
}


void Foam::linearDampingWaveNumber::computeDampingAndFreq
(
    const Foam::scalar REqu
)
{
    REqu_ = REqu;

    p0_ = pInf_ + 2.0*sigma_/REqu;

    scalar Chi = DThermal_/(omega_*sqr(REqu));

    // Complex parameter Phi
    std::complex<scalar> Phi =
        3.0*gamma_
      / (
            1.0 -
            (
                3.0*(gamma_ - 1.0)*i_*Chi
              * (std::sqrt(i_/Chi) * 1.0/std::tanh(std::sqrt(i_/Chi)) - 1.0)
            )
        );

    // Resonant angular and standard frequencies
    omega0Sqr_ =
        p0_/(rhoc_*sqr(REqu)) * (Phi.real() - 2.0*sigma_/(REqu*p0_));

    // Damping factor b
    b_ =
        2.0*mu_ / (rhoc_*sqr(REqu))
      + p0_*Phi.imag() / (2.0*rhoc_*omega_*sqr(REqu))
      + sqr(omega_)*REqu / (2.0*cSound_);
}


Foam::scalar Foam::linearDampingWaveNumber::realPartWithVoidFrac
(
    const Foam::scalar voidFrac
)
{
    return
    (
        sqr(omega_/cSound_)
      + 3.0*voidFrac*sqr(omega_)
      / (sqr(REqu_)*(omega0Sqr_ - sqr(omega_)))
    );
}


Foam::scalar Foam::linearDampingWaveNumber::computeBlakePressure
(
    const Foam::scalar REqu
)
{
    scalar S = 2.0*sigma_/(pInf_*REqu);
    return pInf_*(1.0 + sqrt(4.0/27.0*pow3(S)/(1.0 + S)));
}

// ************************************************************************* //