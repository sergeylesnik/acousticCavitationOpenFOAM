#include "residualNorm.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineRunTimeSelectionTable(residualNorm, dict)
    defineTypeNameAndDebug(residualNorm, 0)
    defineTypeNameAndDebug(scaledL1Norm, 0)
    defineTypeNameAndDebug(normalizedL1Norm, 0)
    defineTypeNameAndDebug(scaledL2Norm, 0)
    defineTypeNameAndDebug(normalizedL2Norm, 0)

    addToRunTimeSelectionTable
    (
        residualNorm,
        scaledL1Norm,
        dict
    );

    addToRunTimeSelectionTable
    (
        residualNorm,
        normalizedL1Norm,
        dict
    );

    addToRunTimeSelectionTable
    (
        residualNorm,
        scaledL2Norm,
        dict
    );

    addToRunTimeSelectionTable
    (
        residualNorm,
        normalizedL2Norm,
        dict
    );
}


Foam::autoPtr<Foam::residualNorm>
Foam::residualNorm::New(const dictionary& dict)
{

    word normType(dict.lookup("residualNorm"));
    dictConstructorTable::iterator cstrIter =
            dictConstructorTablePtr_->find(normType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Foam::residualNorm::New(const dictionary& dict)"
        )   << "Unknown normType "
            << normType
            << ", constructor not in hash table" << nl << nl
            << "    Valid residualNorm types are:" << nl
            << dictConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<residualNorm>(cstrIter()());
}


void Foam::residualNorm::cmpxResidualField
(
   const fvBlockMatrix<vector2> &matrix,
   vector2Field &res
)
{
    // Bring everything to the RHS to obtain -F
    const vector2Field& x = matrix.psi();
    const vector2Field& b = matrix.source();

    // Use given field as scratch space to compute Ax
    matrix.Amul(res, x);

    // residual = -F = b - Ax
    res = b - res;
}


Foam::scalar Foam::scaledL1Norm::residual(const fvBlockMatrix<vector2>& matrix)
{
    const vector2Field& x = matrix.psi();
    label localSize = x.size();
    vector2Field resCmpxField(localSize);
    cmpxResidualField(matrix, resCmpxField);

    // Normalize residual by the number of unknowns.
    label totalSize = returnReduce(localSize, sumOp<label>());

    return gSumMag(resCmpxField)/totalSize;
}


Foam::scalar Foam::normalizedL1Norm::residual(const fvBlockMatrix<vector2> &matrix)
{
    // Residual is the sum of the absolute values of the complex field
    // components (vector2)

    const vector2Field& x = matrix.psi();
    const vector2Field& b = matrix.source();
    vector2Field Ax(x.size());
    matrix.Amul(Ax, x);
    vector2Field resCmpx = b - Ax;

    scalar res = gSumMag(resCmpx);

    // Reference x field consisting of <x> in each entry
    vector2Field xRef(x.size(), gAverage(cmptMag(x)));
    // Reuse resCmpx as scratch space to compute A*xRef
    matrix.Amul(resCmpx, xRef);
    scalar normFactor = gSum(mag(Ax - resCmpx) + mag(b - resCmpx)) + SMALL;

    return res/normFactor;
}


Foam::scalar Foam::scaledL2Norm::residual
(
    const fvBlockMatrix<vector2>& matrix
)
{
    // Residual is the square root of the sum of squared real and imaginary
    // parts.

    label localSize = matrix.psi().size();
    vector2Field resCmpxField(localSize);
    cmpxResidualField(matrix, resCmpxField);

    // Equals (complex)*(complex conjugate)
    scalarField resSqr = cmptSum(cmptMultiply(resCmpxField, resCmpxField));
    scalar res = sqrt(gSum(resSqr));
    // Normalize residual by the number of unknowns.
    label totalSize = returnReduce(localSize, sumOp<label>());

    return res/totalSize;
}

Foam::scalar Foam::normalizedL2Norm::residual
(
    const fvBlockMatrix<vector2>& matrix
)
{
    // Residual is the square root of the sum of squared real and imaginary
    // parts.

    const vector2Field& x = matrix.psi();
    const vector2Field& b = matrix.source();
    vector2Field Ax(x.size());
    matrix.Amul(Ax, x);
    vector2Field resCmpx = b - Ax;

    // Equals (complex)*(complex conjugate)
    scalarField resSqr = cmptSum(cmptMultiply(resCmpx, resCmpx));
    scalar res = sqrt(gSum(resSqr));

    // Normalize residual like in FOAM
    // Reference x field consisting of <x> (averaged x) in each entry
    vector2Field xRef(x.size(), gAverage(x));
    // Reuse resCmpx as scratch space to compute A*xRef
    matrix.Amul(resCmpx, xRef);
    scalar normFactor = gSum(mag(Ax - resCmpx) + mag(b - resCmpx)) + SMALL;

    return res/normFactor;
}
