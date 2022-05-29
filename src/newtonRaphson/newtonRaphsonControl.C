#include "newtonRaphsonControl.H"

Foam::newtonRaphsonControl::newtonRaphsonControl
(
    const dictionary& dict
)
:
    dict_(dict),
    outerLoopResidual_(VGREAT),
    prevTrialResidual_(VGREAT),
    trialResidual_(VGREAT),
    normalizedResidual_(VGREAT),
    errLast_(VGREAT),
    errTry_(VGREAT),
    damping_(2.0),
    prevDamping_(damping_),
    init_(false),
    maxOuterInitIter_(readLabel(dict_.lookup("maxOuterInitIter"))),
    initDamping_(readScalar(dict_.lookup("initDamping"))),
    minDamping_(readScalar(dict_.lookup("minDamping"))),
    maxRecoveryDamping_(readScalar(dict_.lookup("maxRecoveryDamping"))),
    nRaiseRecovery_(readLabel(dict_.lookup("nRaiseRecovery"))),
    nRecoveryStep_(0),
    residualControl_(readScalar(dict_.lookup("residualControl"))),
    maxOuterIter_(readScalar(dict_.lookup("maxOuterIterations"))),
    nInnerIter_(0),
    nOuterIter_(0),
    recoveryStep_(false),
    resNormPtr_(residualNorm::New(dict_)),
    normalizedNorm_()
{}

bool newtonRaphsonControl::outerLoop()
{
    nOuterIter_++;

    scalar relRes = normalizedResidual_;
    bool converged = (residualControl_ > relRes);
    bool completed = (nOuterIter_ > maxOuterIter_);
    if (completed || converged)
    {
        if (converged)
        {
            Info<< nl << "NewtonRaphsonSolver: solution converged in "
                << nOuterIter_ - 1 << " iteration(s)" << nl << endl;
        }
        else  // completed = true
        {
            Info<< nl
                << "NewtonRaphsonSolver: max number of outer iterations = "
                << maxOuterIter_ << " reached. Ending the solver."
                << nl << endl;
        }

        // Reset for the next outer loop
        nOuterIter_ = 0;
        outerLoopResidual_ = VGREAT;
        trialResidual_ = VGREAT;
        normalizedResidual_ = VGREAT;
    }

    return (!converged && !completed);
}

bool newtonRaphsonControl::innerLoop()
{
    nInnerIter_++;

    bool completed = (trialResidual_ < outerLoopResidual_);

    if (completed || recoveryStep_ || init_)
    {
        Info<< "Inner iterations = " << nInnerIter_ - 1
            << ", Damping = " << damping_
            << ", Initial normalized L1 residual = " << normalizedResidual_
            << ", Final user-defined norm residual = " << trialResidual_
            << endl;

        // Reset
        nInnerIter_ = 0;
        damping_ = 2.0;
        init_ = false;

        // Set the current residual to the trial residual from
        // the last inner iteration loop.
        outerLoopResidual_ = trialResidual_;

        if (recoveryStep_)
        {
            recoveryStep_ = false;
            nRecoveryStep_++;
        }
        else
        {
            nRecoveryStep_ = 0;
        }
        return false;
    }

    // Adjust damping
    if (nOuterIter_ <= maxOuterInitIter_)
    {
        // Use initial damping from dictionary at the initializaiton.
        damping_ = initDamping_;
        init_ = true;
    }
    else
    {
        if (nInnerIter_ <= 2)
        {
            prevDamping_ = damping_;
            damping_ /= 2;
        }
        else
        {
            polynomialSearch(damping_, prevDamping_);
        }

        // Set up recovery step when arrived at the min. damping.
        if (damping_ < minDamping_)
        {
            // Gradually raise current damping up to the threshold given by
            // the user.
            scalar currRecovery;

            if (nRecoveryStep_ < nRaiseRecovery_)
            {
                currRecovery =
                    maxRecoveryDamping_
                        / pow(2, nRaiseRecovery_ - nRecoveryStep_);
            }
            else
            {
                currRecovery = maxRecoveryDamping_;
            }

            damping_ = currRecovery;

            recoveryStep_ = true;
            Info<< "Reached minimal damping = " << minDamping_
                << ", Recovery step number = " << nRecoveryStep_
                << ", Current recovery damping = " << damping_  << endl;
        }

    }

    return !completed;
}

void newtonRaphsonControl::polynomialSearch
(
    scalar& damping,
    scalar& prevDamping
)
{
    const scalar d1 = damping;
    const scalar d2 = prevDamping;
    const scalar r0 = outerLoopResidual_;
    const scalar r1 = trialResidual_;
    const scalar r2 = prevTrialResidual_;

    // Polynomial model f
    // f'(0)
    const scalar f1 = (-d2/d1 * (r1 - r0) + d1/d2 * (r2 - r0)) / (d1 - d2);
    // f''(0)
    const scalar f2 =
            2 / (d1*d2 * (d1 - d2)) * (d2 * (r1 - r0) - d1 * (r2 - r0));

    // Save damping for the next iteration.
    prevDamping = damping;

    const scalar sigma0 = 0.1;
    const scalar sigma1 = 0.5;
    if (f2 > 0 && f1 < 0)
    {
        scalar d = -f1 / f2;

        // Safeguarding: don't allow for large variations of damping.
        if (d < sigma0 * d1)
        {
            damping *= sigma0;
        }
        else if (d > sigma1 * d1)
        {
            damping *= sigma1;
        }
        else
        {
            damping = d;
        }
    }
    else
    {
        damping *= sigma1;
    }
}

scalar& newtonRaphsonControl::dampingFactor()
{
    return damping_;
}

void newtonRaphsonControl::updateTrialResidual
(
    const fvBlockMatrix<vector2>& matrix
)
{
    prevTrialResidual_ = trialResidual_;
    trialResidual_ = resNormPtr_->residual(matrix);
}

void newtonRaphsonControl::updateInitialNormalizedResidual
(
    const fvBlockMatrix<vector2>& matrix
)
{
    normalizedResidual_ = normalizedNorm_.residual(matrix);

    // Compute initial user-defined norm residual at the solver start up.
    // In case of solver restart the solver behaviour stays determenistic.
    if (nOuterIter_ == 1)
    {
        outerLoopResidual_ = resNormPtr_->residual(matrix);
    }
}

scalar& newtonRaphsonControl::outerLoopResidual
(
)
{
    return outerLoopResidual_;
}

const dictionary& newtonRaphsonControl::dict()
{
    return dict_;
}
