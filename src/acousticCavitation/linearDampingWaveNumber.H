#ifndef linearDampingWaveNumber_H
#define linearDampingWaveNumber_H

#include <complex>
#include "complex.H"
#include "fvCFD.H"

namespace Foam
{

class linearDampingWaveNumber
{
    // Static Data Members
    
        //- Imaginary unit
        //  Foam::complex doesn't have sqrt() and tanh(); use std::complex<>
        static const std::complex<scalar> i_;
    
    // Private Data
    // Keep both dimensioned and scalar for eventual dimension checks and
    // scalars for speed

        const scalar pi_;
        const dictionary& transportProps_;
        const dictionary& bubbleProps_;
        dictionary& derivedProps_;
        const dimensionedScalar freqDim_;
        const dimensionedScalar pInfDim_;
        const dimensionedScalar REquDim_;
        const dimensionedScalar voidFracDim_;
        const dimensionedScalar rhocDim_;
        const dimensionedScalar cSoundDim_;
        const dimensionedScalar nuDim_;
        const dimensionedScalar sigmaDim_;
        const dimensionedScalar DThermalDim_;
        const dimensionedScalar gammaDim_;
        const dimensionedScalar omegaDim_;
        const dimensionedScalar muDim_;
        scalar kSqrReCavDim_;
        scalar kSqrReNoCavDim_;
        const scalar freq_;
        const scalar pInf_;
        scalar REqu_;
        const scalar voidFrac_;
        const scalar rhoc_;
        const scalar cSound_;
        const scalar nu_;
        const scalar sigma_;
        const scalar DThermal_;
        const scalar gamma_;
        const scalar omega_;
        const scalar mu_;
        const scalar kSqrReNoCav_;
        scalar p0_;
        scalar omega0Sqr_;
        scalar b_;
        scalar kSqrIm_;
        scalar kSqrReCav_;

    // Private Member Functions
    
        //- Compute both real and imaginary parts of the wave number
        void computeComplex(const scalar voidFrac, const scalar REqu);

public:

    // Constructor

        linearDampingWaveNumber
        (
            const dictionary& transportProps,
            const dictionary& bubbleProps,
            dictionary& derivedProps,
            dimensionedScalar REqu = dimensionedScalar("REqu", dimLength, 0.0),
            dimensionedScalar voidFrac =
                dimensionedScalar("voidFrac", dimless, 0.0)
        );
        
    // Member Functions

        //- Compute wave number if both void fraction and eq. radius are set
        void compute();

        //- Compute b and omega0 with given equlibrium radius
        void computeDampingAndFreq(const scalar REqu);
        
        //- Compute only the real part of the wave number
        //  b and omega0 are to be set before by one of the other methods
        scalar realPartWithVoidFrac(const scalar voidFrac);

        //- Compte Blake pressure threshold
        //  Not strictly a part of wave number but all of the variables are already
        //  provided by this class
        scalar computeBlakePressure(const scalar REqu);

    // Access

        inline const dimensionedScalar& rhoLiquidDim() const;

        inline const dimensionedScalar& omegaDim() const;

        inline scalar rhoLiquid() const;

        inline scalar omega() const;

        inline scalar kSqrReCav() const;

        inline scalar kSqrReNoCav() const;

};

} // End namespace Foam


inline const Foam::dimensionedScalar&
Foam::linearDampingWaveNumber::rhoLiquidDim() const
{
    return rhocDim_;
}

inline const Foam::dimensionedScalar&
Foam::linearDampingWaveNumber::omegaDim() const
{
    return omegaDim_;
}

inline Foam::scalar Foam::linearDampingWaveNumber::rhoLiquid() const
{
    return rhoc_;
}

inline Foam::scalar Foam::linearDampingWaveNumber::omega() const
{
    return omega_;
}

inline Foam::scalar Foam::linearDampingWaveNumber::kSqrReCav() const
{
    return kSqrReCav_;
}

inline Foam::scalar Foam::linearDampingWaveNumber::kSqrReNoCav() const
{
    return kSqrReNoCav_;
}

#endif