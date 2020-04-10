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

Description
    Solver for the homogeneous Helmholtz equation using the block-coupled
    matrix class for building of the linear system and a multi-frontal (direct)
    linear solver MUMPS.
    Please refer to the following paper whenever using this code:
    "Lesnik, S., Brenner, G., Ayaz-Bustami, K., & Mettin, R. Modellierung der
    akustischen Kavitation mithilfe des Euler-Lagrange Ansatzes."

Author
    Sergey Lesnik, ITM Clausthal, ITM Clausthal

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "simpleControl.H"

#include "blockMUMPSSolver.H"
#include "blockMUMPSSolverSer.H"
#include "blockMUMPSSolverPar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createFields.H"

    scalar maxResidual = 0;
    scalar convergenceCriterion = 0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting loop\n" << endl;

    // MUMPS needs MPI_Init() to be executed even for a serial run.
    // OF routine doesn't allow this. Use MPI_Init() directly.
    if (!Pstream::parRun())
    {
        MPI_Init(NULL, NULL);
    }

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvScalarMatrix PAcReEqn
        (
            fvm::laplacian(PAcRe)
          + fvm::Sp(kSqrRe, PAcRe)
          - kSqrIm*PAcIm
        );

        fvScalarMatrix PAcImEqn
        (
            fvm::laplacian(PAcIm)
          + fvm::Sp(kSqrRe, PAcIm)
          + kSqrIm*PAcRe
        );

        // Prepare block matrix
        fvBlockMatrix<vector2> helmholtzEqn(blockT);

        // Insert equations into block matrix
        helmholtzEqn.insertEquation(0, PAcReEqn);
        helmholtzEqn.insertEquation(1, PAcImEqn);

        // Add off-diagonal coupling terms
        scalarField coupling(mesh.nCells(), kSqrIm.value());
        helmholtzEqn.insertEquationCoupling(0, 1, -coupling);
        helmholtzEqn.insertEquationCoupling(1, 0, coupling);

        // Update source coupling: coupling terms eliminated from source
        helmholtzEqn.updateSourceCoupling();

        // Initialize MUMPS using OF's factory method
        // Analyzation and factorization phases are also done here
        autoPtr<blockMUMPSSolver> directSolver =
                blockMUMPSSolver::New(helmholtzEqn, mesh);

        directSolver->getRhs();
        directSolver->dumpAccToMUMPSDict();

        // Runs always at least ones (e.g. for nNonOrthogonalCorrectors 0;)
        while (simple.correctNonOrthogonal())
        {
            Info<< "Non-orthogonal correction = "
                << simple.corrNonOrtho() - 1 << nl << endl;

            PAcRe.storePrevIter();
            PAcIm.storePrevIter();

            fvScalarMatrix PAcReEqnCorr
            (
                fvm::laplacian(PAcRe)
              + fvm::Sp(kSqrRe, PAcRe)
              - kSqrIm*PAcIm
            );

            fvScalarMatrix PAcImEqnCorr
            (
                fvm::laplacian(PAcIm)
              + fvm::Sp(kSqrRe,PAcIm)
              + kSqrIm*PAcRe
            );

            // Build a new fvBlockMatrix with the corrected RHS.
            // Work around for the problem: insertion into the first
            // matrix (helmholtzEqn) object leads to doubling of the RHS.
            fvBlockMatrix<vector2> helmholtzEqnCorr(blockT);
            helmholtzEqnCorr.insertEquation(0, PAcReEqnCorr);
            helmholtzEqnCorr.insertEquation(1, PAcImEqnCorr);
            helmholtzEqnCorr.insertEquationCoupling(0, 1, -coupling);
            helmholtzEqnCorr.insertEquationCoupling(1, 0, coupling);
            helmholtzEqnCorr.updateSourceCoupling();

            // Solve with the RHS from the corrected matrix.
            // Note that the solution will be saved in the matrix the solver
            // was initialized with (helmholtzEqn).
            directSolver->solveWithRhs(helmholtzEqnCorr);

            helmholtzEqnCorr.retrieveSolution(0, PAcRe.internalField());
            helmholtzEqnCorr.retrieveSolution(1, PAcIm.internalField());

            PAcRe.correctBoundaryConditions();
            PAcIm.correctBoundaryConditions();

            // Relax using solution either from the previous time step which
            // was written by simple.loop() or from .prevIter() if
            // storePrevIter() has been called. Since no solution is written
            // while the non-orthogonal correction loop, call storePrevIter()
            // in the beginning of the loop.
            PAcRe.relax();
            PAcIm.relax();

            PAc = sqrt(pow(PAcRe,2) + pow(PAcIm,2));
            scalar maxPAc = gMax(PAc);

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << "  Max(PAc) = " << maxPAc << " Pa"
                << nl << endl;

            // Check convergence
            mesh.solutionDict().subDict("blockSolver").readIfPresent
            (
                "convergence",
                convergenceCriterion
            );

            maxResidual = cmptMax(directSolver->initialResidual());

            if (maxResidual < convergenceCriterion)
            {
                Info<< "reached convergence criterion: "
                    << convergenceCriterion << endl;

                // Ends the simple loop only. Use 'break' to end the
                // non-orthogonal correction loop.
                runTime.writeAndEnd();
                break;
            } // End non-othogonal corrections loop

        } // End SIMPLE loop

    } // End main


    if (!Pstream::parRun())
    {
        MPI_Finalize();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
