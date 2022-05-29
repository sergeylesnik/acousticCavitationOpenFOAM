#ifndef waveNumTableInterpolator_H
#define waveNumTableInterpolator_H

#include "fvCFD.H"
#include "interpolation2DTable.H"
#include "AcCavitationCloud.H"
#include "stepSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class waveNumCalc
:
    protected stepSmoother
{

protected:

    // Protected Data
    const fvMesh& mesh_;
    const dictionary& dict_;

public:

    // Constructor

        waveNumCalc
        (
            const fvMesh& mesh,
            const dictionary& dict
        )
        :
            stepSmoother
            (
                dict.lookupOrDefault("cavitationInceptionSmoothing", 0.03)
            ),
            mesh_(mesh),
            dict_(dict)
        {}

    // Dectructor

        virtual ~waveNumCalc(){}

    // Member Functions

        virtual void compute
        (
            volScalarField& kSqrIm,
            volScalarField& kSqrRe,
            const volScalarField& P
        ) = 0;
};


class waveNumCloudInterpolator
:
public waveNumCalc
{

protected:

    AcCavitationCloud<acCavitationParcel>& bubbleCloud_;
    const dimensionedScalar rho_;
    linearDampingWaveNumber& ldwn_;
    const volScalarField& voidFrac_;

public:

    // Constructors
    waveNumCloudInterpolator
    (
        fvMesh& mesh,
        AcCavitationCloud<acCavitationParcel>& cloud,
        const dimensionedScalar rho,
        linearDampingWaveNumber& ldwn,
        const volScalarField& voidFrac
    )
    :
        waveNumCalc(mesh, cloud.constProps().dict()),
        bubbleCloud_(cloud),
        rho_(rho),
        ldwn_(ldwn),
        voidFrac_(voidFrac)
    {}

    // Member functions
    virtual void compute
    (
        volScalarField& kSqrIm,
        volScalarField& kSqrRe,
        const volScalarField& P
    )
    {
        bubbleCloud_.computeImagWaveNumber
        (
            kSqrIm.dimensionedInternalField(),
            kSqrRe.dimensionedInternalField(),
            P.dimensionedInternalField(),
            voidFrac_.dimensionedInternalField(),
            ldwn_,
            *this
        );
    }

};


class waveNumTableInterpolator
:
public waveNumCalc
{

protected:

    const interpolation2DTable<scalar>& nonLinDampingTable_;
    const dimensionedScalar REqu_;
    const dimensionedScalar kSqrRePure_;
    const dimensionedScalar kSqrReCav_;
    const dimensionedScalar sigma_;
    const dimensionedScalar pInf_;
    dimensionedScalar PBlake_;
    dimensionedScalar lowerP_;
    dimensionedScalar upperP_;

public:

    // Constructors

        waveNumTableInterpolator
        (
            const interpolation2DTable<scalar>& nonLinDampingTable,
            const fvMesh& mesh,
            const dictionary& dict
        )
        :
            waveNumCalc(mesh, dict),
            nonLinDampingTable_(nonLinDampingTable),
            REqu_(dict_.lookup("REqu")),
            kSqrRePure_(dict_.lookup("kSqrRePure")),
            kSqrReCav_(dict_.lookup("kSqrReCav")),
            sigma_(dict_.lookup("sigma")),
            pInf_(dict_.lookup("pInf")),
            PBlake_(pInf_),
            lowerP_(pInf_),
            upperP_(pInf_)
        {
            dimensionedScalar S = 2.0 * sigma_ / (pInf_ * REqu_);
            PBlake_ = pInf_ * (1.0 + sqrt(4.0/27.0 * pow(S, 3) / (1.0 + S)));
            lowerP_ = (1.0 - spread_) * PBlake_;
            upperP_ = (1.0 + spread_) * PBlake_;
        }

    // Member Functions

        void clampAndSmooth
        (
            volScalarField& kSqrIm,
            volScalarField& kSqrRe,
            const volScalarField& P
        )
        {
            forAll(mesh_.C(), cellI)
            {
                scalar nuclCoeff =
                    smoothStep
                    (
                        lowerP_.value(), upperP_.value(), P[cellI]
                    );
                kSqrIm[cellI] *= nuclCoeff;
                kSqrRe[cellI] =
                    (1.0 - nuclCoeff) * kSqrRePure_.value()
                  + nuclCoeff * kSqrReCav_.value();
            }
        }

        void saveDerivedParamsToDict(IOdictionary& dict)
        {
            dict.add("PAcBlake", PBlake_);
            dict.add("lowerPAcSmoothingLimit", lowerP_);
            dict.add("upperPAcSmoothingLimit", upperP_);
        }

        virtual void compute
        (
            volScalarField& kSqrIm,
            volScalarField& kSqrRe,
            const volScalarField& P
        )
        {
            forAll(mesh_.C(), cellI)
            {
                kSqrIm[cellI] = nonLinDampingTable_(REqu_.value(), P[cellI]);
            }
            clampAndSmooth(kSqrIm, kSqrRe, P);
        }

};

// Introduces inhomogeneous void fraction.
class waveNumTableVoidFracInterpolator
:
public waveNumTableInterpolator
{
private:

    const volScalarField& voidFrac_;
    const scalar voidFrac0_;
    const dimensionedScalar omega_;
    const dimensionedScalar cSound_;

public:

    // Constructor

        waveNumTableVoidFracInterpolator
        (
            const interpolation2DTable<scalar>& nonLinDampingTable,
            const fvMesh& mesh,
            const volScalarField& voidFrac,
            const scalar voidFrac0,
            const dictionary& dict
        )
        :
            waveNumTableInterpolator(nonLinDampingTable, mesh, dict),
            voidFrac_(voidFrac),
            voidFrac0_(voidFrac0),
            omega_
            (
                dimensionedScalar
                (
                    mesh_.lookupObject<dictionary>
                    (
                        "derivedProperties"
                    ).lookup("omega")
                )
            ),
            cSound_((dimensionedScalar(dict_.lookup("cSound"))))
        {}

    // Member Functions
    
        virtual void compute
        (
            volScalarField& kSqrIm,
            volScalarField& kSqrRe,
            const volScalarField& P
        )
        {
            forAll(mesh_.C(), cellI)
            {
                kSqrIm[cellI] = nonLinDampingTable_(REqu_.value(), P[cellI]);

                // Introduce inhomogeneous void fraction
                kSqrIm[cellI] *= voidFrac_[cellI] / voidFrac0_;
                dimensionedScalar subterm = pow(omega_, 2) / pow(cSound_, 2);
                kSqrRe[cellI] =
                    (kSqrRe[cellI] - subterm.value())
                  * voidFrac_[cellI] / voidFrac0_ + subterm.value();
            }

            clampAndSmooth(kSqrIm, kSqrRe, P);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // waveNumTableInterpolator_H

// ************************************************************************* //
