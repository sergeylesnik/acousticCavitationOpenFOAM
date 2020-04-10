/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blockMUMPSSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockMUMPSSolver, 0);
    defineRunTimeSelectionTable(blockMUMPSSolver, runType)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMUMPSSolver::blockMUMPSSolver
(
    fvBlockMatrix<vector2>& matrix,
    const fvMesh& mesh
)
:
    matrix_(matrix),
    mesh_(mesh),
    nCells_(mesh_.nCells()),
    nInternalFaces_(matrix_.lduAddr().lowerAddr().size()),
    irn_(1),
    jcn_(1),
    amv_(1),
    rhs_(1),
    myid_(0),
    matrixI_(0),
    initialNormResidual_(pTraits<vector2>::zero),
    finalNormResidual_(pTraits<vector2>::zero),
    MUMPSdict_(mesh_.lookupObject<IOdictionary>("MUMPSSettings")),
    dumpCompleteLinSys_
    (
        MUMPSdict_.lookupOrDefault<Switch>("dumpCompleteLinSys", 0)
    ),
    dumpFilePrefix_
    (
        MUMPSdict_.lookupOrDefault<word>("dumpFilePrefix", "")
    ),
    printResiduals_
    (
        MUMPSdict_.lookupOrDefault<Switch>("printResiduals", 1)
    )
{
    // Initialize a MUMPS instance. Use MPI_COMM_WORLD
    mumps_.job = JOB_INIT;
    mumps_.par = 1;  // non-working (0) or working (1) master
    mumps_.sym = 0;  // non-symmetric (0) or symmetric (1) matrix
    mumps_.comm_fortran = USE_COMM_WORLD;

    // Initialize MUMPS with double precision arithmetics (d)
    dmumps_c(&mumps_);

    // Get and set MUMPS settings from IOdictionary
    HashTable<label, label> icntlTable(MUMPSdict_.lookup("ICNTL"));
    for
    (
        HashTable<label, label>::iterator icntlIter = icntlTable.begin();
        icntlIter != icntlTable.end();
        ++icntlIter
    )
    {
        label settingI = icntlIter.key();
        mumps_.ICNTL(settingI) = icntlTable[settingI];
    }

    bool dumpByMUMPS
    (
        MUMPSdict_.lookupOrDefault<Switch>("dumpLinSysByMUMPS", 0)
    );
    if (dumpByMUMPS)
    {
        const word& fileName(MUMPSdict_.lookup("dumpByMUMPSFileName"));
        strcpy(mumps_.write_problem, fileName.c_str());
    }

    // Format of RHS: 0 (dense RHS)
    mumps_.ICNTL(20) = 0;

    // Centralized or distributed solution. Set to centralized (0) because
    // distribution is done to MUMPS mapping which is of no use.
    mumps_.ICNTL(21) = 0;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blockMUMPSSolver>
Foam::blockMUMPSSolver::New
(
    fvBlockMatrix<vector2>& matrix,
    const fvMesh& mesh
)
{
    runTypeConstructorTable::iterator cstrIter;

    if (Pstream::parRun())
    {
        cstrIter =
            runTypeConstructorTablePtr_->find("blockMUMPSSolverPar");
    }
    else
    {
        cstrIter =
            runTypeConstructorTablePtr_->find("blockMUMPSSolverSer");
    }

    return autoPtr<blockMUMPSSolver>(cstrIter()(matrix, mesh));
}

void Foam::blockMUMPSSolver::solve()
{
    solveWithRhs(matrix_);
}

void Foam::blockMUMPSSolver::solveWithRhs
(
    const fvBlockMatrix<Foam::vector2>& matrix
)
{
    const Field<vector2>& x = matrix.psi().field();
    const Field<vector2>& b = matrix.source();
    initialNormResidual_ = computeResidual(x, b, matrix);

    getRhs(matrix);

    solveCore();

    passSol();

    finalNormResidual_ = computeResidual(x, b, matrix_);

    if (printResiduals_)
    {
        printResiduals();
    }
}

Foam::vector2 Foam::blockMUMPSSolver::computeResidual
(
    const Field<vector2>& x,
    const Field<vector2>& b,
    const fvBlockMatrix<vector2>& matrix
)
{
    Field<vector2> Ax(rhs_[myid_].size());

    // Residual is the sum of the absolute values of the residual
    // vector components
    matrix.Amul(Ax, x);
    Field<vector2> resVec = b - Ax;
    vector2 res = gSumCmptMag(resVec);

    // Normalisation factor
    Field<vector2> xRef(x.size(), gAverage(x));
    Field<vector2> bRef;
    matrix.Amul(bRef, xRef);
    vector2 normFactor =
        gSumCmptMag(Ax - bRef) + gSumCmptMag(b - bRef) + vector2(SMALL);

    return cmptDivide(res, normFactor);
}

void Foam::blockMUMPSSolver::printResiduals()
{
    Info<< "MUMPS:  Solving for (PAcRe PAcIm)"
        << ", Initial residual = " << initialNormResidual_
        << ", Final residual = " << finalNormResidual_
        << endl;
}

void Foam::blockMUMPSSolver::analyzeAndFactorizeMUMPS()
{
    // job=4 combines job=1 (analyze) job=2 (factorize)
    mumps_.job = 4;
    dmumps_c(&mumps_);
}

void Foam::blockMUMPSSolver::solveMUMPS()
{
    // job=3 computes the solution
    mumps_.job = 3;
    dmumps_c(&mumps_);
}

void Foam::blockMUMPSSolver::assembleDiag()
{
    label j;
    for (label i=0; i<nCells_; i++)
    {
        j = diagInd(i);

        // Note that the indices of the C array start with 0 but
        // the indices of the matrix start with 1 (acc. to MUMPS, Fortran)
        addMatrixEntry(j+1, j+1, matrix_.diag().asSquare()[i][0]);
        addMatrixEntry(j+1, j+2, matrix_.diag().asSquare()[i][1]);
        addMatrixEntry(j+2, j+1, matrix_.diag().asSquare()[i][2]);
        addMatrixEntry(j+2, j+2, matrix_.diag().asSquare()[i][3]);
    }
}

void Foam::blockMUMPSSolver::assembleOffDiag()
{
    label ownerI;
    label neighbourI;

    for (label faceI=0; faceI<nInternalFaces_; faceI++)
    {
        setOwnerNeighbInd(ownerI, neighbourI, faceI);

        // The coefficients of the off-diagonal terms are saved asLinear
        // (meaning as vectors), since only laplace operator discretisation
        // is treated here and these terms are not "block-coupled".
        // Thus, the entries of upper() and lower() are the diagonal terms
        // of the mini-blocks from the block-coupled system

        // Note that the indices of the C array start with 0
        // but the indices of the matrix start with 1 (acc. to MUMPS).
        // Consider only the upper matrix since it's equal to the lower.
        addMatrixEntry
        (
            ownerI+1,
            neighbourI+1,
            matrix_.upper().asLinear()[faceI][0]
        );

        addMatrixEntry
        (
            ownerI+2,
            neighbourI+2,
            matrix_.upper().asLinear()[faceI][1]
        );
        addMatrixEntry
        (
            neighbourI+1,
            ownerI+1,
            matrix_.upper().asLinear()[faceI][0]
        );
        addMatrixEntry
        (
            neighbourI+2,
            ownerI+2,
            matrix_.upper().asLinear()[faceI][1]
        );
    }
}

void Foam::blockMUMPSSolver::addMatrixEntry
(
    label rowI, label colI, scalar aCoeff
)
{

    irn_[myid_][matrixI_] = rowI;
    jcn_[myid_][matrixI_] = colI;
    amv_[myid_][matrixI_] = aCoeff;

    matrixI_++;
}

void Foam::blockMUMPSSolver::getRhs()
{
    getRhs(matrix_);
}

void Foam::blockMUMPSSolver::getRhs(const fvBlockMatrix<vector2>& matrix)
{
    for (label i=0; i<nCells_; i++)
    {
        rhs_[myid_][i*2] = matrix.source()[i][0];
        rhs_[myid_][i*2+1] = matrix.source()[i][1];
    }
}

void Foam::blockMUMPSSolver::passSol()
{
    vector2Field& x = matrix_.psi();
    forAll(x, cellI)
    {
        x[cellI][0] = rhs_[myid_][cellI*2];
        x[cellI][1] = rhs_[myid_][cellI*2+1];
    }
}

void Foam::blockMUMPSSolver::dumpMasterLinearSystem
(
    fileName filePrefix
)
{
    if (Pstream::parRun())
    {
        Pstream::gatherList(irn_);
        Pstream::gatherList(jcn_);
        Pstream::gatherList(amv_);
    }

    if (myid_ == 0)
    {
        word fileName = filePrefix + "Matrix_Time" + mesh_.time().timeName();

        if (debug)
        {
            Info<< "Dumping complete linear system by master proc to "
                << fileName << endl;
        }

        OFstream osMatrix(fileName);

        // Header
        osMatrix << "procNo localI globRow globCol coeff" << endl;

        forAll(irn_, procI)
        {
            forAll(irn_[procI], i)
            {
                osMatrix.precision(12);
                osMatrix << procI << " " << i
                   << " " << irn_[procI][i]
                   << " " << jcn_[procI][i]
                   << " " << amv_[procI][i] << endl;
            }
        }
    }

    dumpCompleteRhs(filePrefix);
}
