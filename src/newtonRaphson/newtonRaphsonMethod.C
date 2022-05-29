
#include "newtonRaphsonMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "residualNorm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newtonRaphsonMethod::newtonRaphsonMethod
(
    newtonRaphsonControl& ctrl,
    autoPtr<blockMUMPSSolver>& directSolver,
    fvBlockMatrix<vector2>& jacobian,
    volScalarField& PAcRe,
    volScalarField& PAcIm,
    volScalarField& kSqrRe,
    volScalarField& kSqrIm,
    waveNumCalc& kSqrInterp
)
:
    ctrl_(ctrl),
    linearSolver_(directSolver),
    mesh_(linearSolver_->mesh()),
    jacobianEqn_(jacobian),
    PAcRe_(PAcRe),
    PAcIm_(PAcIm),
    kSqrRe_(kSqrRe),
    kSqrIm_(kSqrIm),
    kSqrInterp_(kSqrInterp),
    blockPAc_
    (
        IOobject
        (
            "blockF0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector2("zero", dimless, vector2::zero)
    )
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::newtonRaphsonMethod::solve
(
    volScalarField& dPAcRe,
    volScalarField& dPAcIm
)
{
    // Save fields to be able to reset while inner iterations.
    PAcRe_.storePrevIter();
    PAcIm_.storePrevIter();

    // Find Newton direction. Jacobian is already factorized (in the
    // constructor of blockMUMPSSolver). Exchange the right hand side and
    // solve.
    fvBlockMatrix<vector2> helmholtzEqn(blockPAc_);
    buildEqnsAndInsert(helmholtzEqn);

    ctrl_.updateInitialNormalizedResidual(helmholtzEqn);

    Field<vector2>& b = helmholtzEqn.source();
    Field<vector2>& x = helmholtzEqn.psi();
    Field<vector2> Ax(x.size());
    helmholtzEqn.Amul(Ax, x);
    b -= Ax;
    linearSolver_->getRhs(helmholtzEqn);
    linearSolver_->dumpAccToMUMPSDict();
    linearSolver_->solveWithRhs(helmholtzEqn);

    jacobianEqn_.retrieveSolution(0, dPAcRe.internalField());
    jacobianEqn_.retrieveSolution(1, dPAcIm.internalField());

    while
    (
        ctrl_.innerLoop()
    )
    {
        // Update solution. Use P fields from the outer loop (prevIter()).
        PAcRe_ = PAcRe_.prevIter() + ctrl_.dampingFactor()*dPAcRe;
        PAcIm_ = PAcIm_.prevIter() + ctrl_.dampingFactor()*dPAcIm;
        PAcRe_.correctBoundaryConditions();
        PAcIm_.correctBoundaryConditions();

        // Evaluate residual
        // Update wave number. Needed only for residual computation.
        volScalarField PAc(sqrt(pow(PAcRe_,2) + pow(PAcIm_,2)));
        kSqrInterp_.compute(kSqrIm_, kSqrRe_, PAc);
        // Values in blockPAc_ will be replaced acc. to discretization.
        // Thus, no reset of the field needed.
        fvBlockMatrix<vector2> hhEqnInner(blockPAc_);
        buildEqnsAndInsert(hhEqnInner);

        ctrl_.updateTrialResidual(hhEqnInner);
    }

}


void Foam::newtonRaphsonMethod::buildEqnsAndInsert
(
    fvBlockMatrix<vector2>& hhEqn
)
{
    fvScalarMatrix PAcReEqn
    (
        fvm::laplacian(PAcRe_)
        + fvm::Sp(kSqrRe_, PAcRe_)
        - kSqrIm_*PAcIm_
    );

    fvScalarMatrix PAcImEqn
    (
        fvm::laplacian(PAcIm_)
        + fvm::Sp(kSqrRe_, PAcIm_)
        + kSqrIm_*PAcRe_
    );

    // Insert equations into block matrix
    hhEqn.insertEquation(0, PAcReEqn);
    hhEqn.insertEquation(1, PAcImEqn);
}


Foam::scalar
newtonRaphsonMethod::scaledL2NormError(fvBlockMatrix<vector2> &matrix)
{

    // Build -F(x_k+1)
    const Field<vector2>& x = matrix.psi();
    Field<vector2>& b = matrix.source();
    Field<vector2> Ax(x.size());
    matrix.Amul(Ax, x);
    b -= Ax;

    // F'(x_k)e = -F(x_k+1)
    linearSolver_->solveWithRhs(matrix);
    const Field<vector2>& e = jacobianEqn_.psi();

    // (complex Field)*(complex conjugate Field)
    Field<scalar> eSqr = cmptSum(cmptMultiply(e, e));
    label totalSize = mesh_.globalData().nTotalCells();
    scalar normErr = sqrt(gSum(eSqr) / totalSize);

    return normErr;
}


Foam::scalar Foam::newtonRaphsonMethod::weightedL2NormError
(
    const Field<vector2>& e,
    const Field<vector2>& x
)
{
    // Weighted error. Take care about division by 0.
    vector2 vecSmall = vector2::one * SMALL;
    Field<vector2> wE = cmptDivide(e, x + vecSmall);

//    Info<< "e = " << e << nl << "x = " << x << endl;

    // Reuse wE
    wE = cmptMultiply(wE, wE);
    label nCells = x.size();

    return sqrt(gSum(cmptSum(wE)) / (2*nCells));
}
