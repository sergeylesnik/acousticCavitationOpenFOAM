#ifndef numericalDifferential_H
#define numericalDifferential_H

#include "fvCFD.H"
#include "interpolation2DTable.H"
#include "waveNumTableInterpolator.H"
#include "complexScalarField.H"
#include "IOmanip.H"

class numericalDifferential
{
private:

    // Private Data
    const fvMesh& mesh_;
    const volScalarField& P_;
    const dimensionedScalar pInf_;
    const dimensionedScalar disturbancePressure_;
    volScalarField dP_;
    dimensionedScalar PSmall_;

public:

    // Constructors

        numericalDifferential
        (
            const volScalarField& P,
            const fvMesh& mesh,
            const dictionary& dict
        )
        :
            mesh_(mesh),
            P_(P),
            pInf_(dict.lookup("pInf")),
            disturbancePressure_
            (
                dict.lookupOrDefault
                (
                    "dP",
                    dimensionedScalar("dP", dimPressure, 0.1)
                )
            ),
            dP_(P_),
            PSmall_("PSmall", dimPressure, VSMALL)
        {
            dP_ = disturbancePressure_;
        }

    // Member Functions

        volScalarField centralDiffWaveNumber
        (
            volScalarField& kSqrIm,
            volScalarField& kSqrRe,
            waveNumTableVoidFracInterpolator& kSqrInterp
        )
        {
            volScalarField PDistPlus(P_ + dP_);
            volScalarField PDistMinus(P_ - dP_);
            volScalarField kSqrImPlus(kSqrIm);
            volScalarField kSqrImMinus(kSqrIm);

            kSqrInterp.compute(kSqrImPlus, kSqrRe, PDistPlus);
            kSqrInterp.compute(kSqrImMinus, kSqrRe, PDistMinus);

            return (kSqrImPlus - kSqrImMinus) / (2*dP_ + PSmall_);
        }

        void backwardDiffWaveNumber
        (
            volScalarField& kSqrRe,
            volScalarField& kSqrIm,
            const complexScalarField& PCmpx,
            waveNumTableVoidFracInterpolator& kSqrInterp,
            volScalarField& dKidPr,
            volScalarField& dKrdPr,
            volScalarField& dKidPi,
            volScalarField& dKrdPi
        )
        {
            volScalarField PReDist(PCmpx.Re() + dP_);
            volScalarField PImDist(PCmpx.Im() + dP_);
            volScalarField kSqrReDist(kSqrRe);
            volScalarField kSqrImDist(kSqrIm);

            // disturb PRe
            volScalarField PDist = sqrt(pow(PReDist,2) + pow(PCmpx.Im(),2));
            kSqrInterp.compute(kSqrImDist, kSqrReDist, PDist);
            dKidPr = (kSqrImDist - kSqrIm) / (dP_);
            dKrdPr = (kSqrReDist - kSqrRe) / (dP_);

            // disturb PIm
            PDist = sqrt(pow(PCmpx.Re(),2) + pow(PImDist,2));
            kSqrInterp.compute(kSqrImDist, kSqrReDist, PDist);
            dKidPi = (kSqrImDist - kSqrIm) / (dP_);
            dKrdPi = (kSqrReDist - kSqrRe) / (dP_);
        }

        // Uses either forward or backward differencing scheme depending on the
        // sign of the acoustic pressure value. Differencing step dP is variable
        // depending on the value of P. The step value is kept sufficiently
        // small to approximate differential and sufficiently high to prevent
        // numerical cancelation errors.
        // Reference: Dennis, Schnabel 1996, pp.98
        void firstOrderDiffWaveNumber
        (
            volScalarField& kSqrRe,
            volScalarField& kSqrIm,
            const complexScalarField& PCmpx,
            waveNumCalc& kSqrInterp,
            volScalarField& dKidPr,
            volScalarField& dKrdPr,
            volScalarField& dKidPi,
            volScalarField& dKrdPi
        )
        {
            // Work with the internal fields only
            scalarField& dKidPrI = dKidPr.internalField();
            scalarField& dKrdPrI = dKrdPr.internalField();
            scalarField& dKidPiI = dKidPi.internalField();
            scalarField& dKrdPiI = dKrdPi.internalField();
            scalarField& PCmpxReI = PCmpx.Re().internalField();
            scalarField& PCmpxImI = PCmpx.Im().internalField();
            scalarField& kSqrReI = kSqrRe.internalField();
            scalarField& kSqrImI = kSqrIm.internalField();
            
            // Define steps for each element of the fields.
            // max(field, scalar) is overloaded to deliver a field with
            // max(field[i], scalar) for each element of the field.
            scalarField dPRe =
                    ROOTSMALL * max(cmptMag(PCmpxReI), pInf_.value())
                  * sign(PCmpxReI);
            // Info<< "dPRe = " << dPRe << endl;
            scalarField dPIm =
                    ROOTSMALL * max(cmptMag(PCmpxImI), pInf_.value())
                  * sign(PCmpxImI);

            // A trick which assures that (P + dP) field consist of exact
            // computer numbers. Taken from Dennis, Schnabel 1996, p.99.
            dPRe += PCmpxReI;
            dPRe -= PCmpxReI;
            dPIm += PCmpxImI;
            dPIm -= PCmpxImI;

            volScalarField kSqrReDist(kSqrRe);
            volScalarField kSqrImDist(kSqrIm);
            scalarField& kSqrReDistI = kSqrReDist.internalField(); 
            scalarField& kSqrImDistI = kSqrImDist.internalField(); 
            volScalarField PDist(PCmpx.Re());
            scalarField& PDistI = PDist.internalField();

            // disturb PRe
            PDistI = PCmpxReI + dPRe;
            PDistI = sqrt(sqr(PDistI) + sqr(PCmpxImI));
            kSqrInterp.compute(kSqrImDist, kSqrReDist, PDist);
            dKidPrI = (kSqrImDistI - kSqrImI) / dPRe;
            dKrdPrI = (kSqrReDistI - kSqrReI) / dPRe;

            // disturb PIm
            PDistI = PCmpxImI + dPIm;
            PDistI = sqrt(sqr(PCmpxReI) + sqr(PDistI));
            kSqrInterp.compute(kSqrImDist, kSqrReDist, PDist);
            dKidPiI = (kSqrImDistI - kSqrImI) / dPIm;
            dKrdPiI = (kSqrReDistI - kSqrReI) / dPIm;
            
            // Correct boundary values since only internal field had been 
            // computed here. Especially processor boundaries might be
            // initialized to garbage values.
            dKidPr.correctBoundaryConditions();
            dKrdPr.correctBoundaryConditions();
            dKidPi.correctBoundaryConditions();
            dKrdPi.correctBoundaryConditions();
        }

};

#endif // numericalDifferential_H

