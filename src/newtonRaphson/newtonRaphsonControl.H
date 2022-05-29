#ifndef newtonRaphsonControl_H
#define newtonRaphsonControl_H

#include "fvCFD.H"
#include "blockMUMPSSolver.H"
#include "residualNorm.H"

namespace Foam
{

class newtonRaphsonControl
{
    const dictionary& dict_;

    // Residuals and errors
    scalar outerLoopResidual_;
    scalar prevTrialResidual_;
    scalar trialResidual_;
    scalar normalizedResidual_;
    scalar errLast_;
    scalar errTry_;

    // Damping coefficients and damping controls
    scalar damping_;
    scalar prevDamping_;
    bool init_;
    label maxOuterInitIter_;
    scalar initDamping_;
    scalar minDamping_;
    scalar maxRecoveryDamping_;
    label nRaiseRecovery_;
    label nRecoveryStep_;

    // Convergence control
    scalar residualControl_;
    label maxOuterIter_;
    label nInnerIter_;
    label nOuterIter_;
    bool recoveryStep_;

    // Residual norm
    autoPtr<residualNorm> resNormPtr_;
    normalizedL1Norm normalizedNorm_;

public:

    // Constructors

        newtonRaphsonControl(const dictionary& dict);

    // Member functions

        bool outerLoop();

        bool innerLoop();

        // A backtracking algorithm. Searches for the best suitable damping.
        // Three-point parabolic line search uses 3 last residual/error
        // values to approximate the damping factor which corresponds to the
        // minimum of the modeled residual/error.
        // Reference: Kelley 1995 p.143
        void polynomialSearch(scalar& damping, scalar& prevDamping);

    // Access

        scalar& dampingFactor();

        void updateTrialResidual(const fvBlockMatrix<vector2>& matrix);

        void updateInitialNormalizedResidual(const fvBlockMatrix<vector2>& matrix);

        scalar& outerLoopResidual();

        const dictionary& dict();
};

} // End namespace Foam

#endif // newtonRaphsonControl_H
