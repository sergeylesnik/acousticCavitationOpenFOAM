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

Application
    acousticCavitationCloudFoam

Description
    The solver includes the solution of the Helmholtz equation, discrete
    cavitation bubbles and URANS modeling of the surrounding liquid. The effect
    of the oscillating bubbles on the acoustics is achieved by introducing the
    attenuation of the acoustic field due to losses during the bubble
    oscillations. The latter are computed using 2D interpolation tables obtained
    from a bubble radial dynamics solver in a pre-processing step. A similar
    approach is used to depict the effect of the acoustic waves on the bubble
    motion (primary Bjerknes force). The coupling between the bubbles and the
    liquid is treated with standard OpenFOAM routines.

Author
    Sergey Lesnik, 2022

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvBlockMatrix.H"

// Multiphase
#include "AcCavitationCloud.H"
#include "acCavitationParcel.H"

// Direct solver
#include "blockMUMPSSolver.H"
#include "blockMUMPSSolverSer.H"
#include "blockMUMPSSolverPar.H"

// Load balancing
#include "loadBalanceFvMesh.H"

// Others
#include "interpolation2DTable.H"
#include "numericalDifferential.H"
#include "waveNumTableInterpolator.H"
#include "newtonRaphsonMethod.H"
#include "absAndPhaseOfVectorField.H"
#include "loadBalanceChecker.H"
#include "linearDampingWaveNumber.H"

#include "bugFixedPimpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    bugFixedPimpleControl pimple(mesh);

#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "createControls.H"


    newtonRaphsonControl loopControl(ctrlDict);

    loadBalanceChecker loadChecker(mesh, bubbleCloud);

    voidFrac = bubbleCloud.theta();

    waveNumCloudInterpolator kSqrCloudInterp
    (
        mesh,
        bubbleCloud,
        rho,
        kSqrLinDamp,
        voidFrac
    );

    if (!Pstream::parRun())
    {
        // MUMPS needs MPI_Init() to be executed even for a serial run.
        // OF routine doesn't allow this. Use MPI_Init() directly.
        MPI_Init(NULL, NULL);
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Load balancing
        bool meshChanged = false;
        if (loadChecker.isImba())
        {
            meshChanged = mesh.update();
            reduce(meshChanged, orOp<bool>());
        }

        //// Lagrangian
        //
        // Initial injection accomplished at object construction
        bubbleCloud.evolve();
        voidFrac = bubbleCloud.theta();
        voidFrac.correctBoundaryConditions();


        // blockT is used as dummy: sets blockJac to zero every iteration.
        // Resetting directly would be less efficient (?).
        // Used only for residual computation in the nonLin solver.
        // Creation of blockT has to happen here (not in createFields) because
        // loadBalanceFvMesh doesn't redistribute volVectorNFields.
        volVector2Field blockT
        (
            IOobject
            (
                "blockT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector2("zero", dimless, vector2::zero)
        );

        //// Newton-Raphson
        while (loopControl.outerLoop())
        {
            volVector2Field blockJac(blockT);

            kSqrCloudInterp.compute(kSqrIm, kSqrRe, PAc);

            complexScalarField PCmpx(PAcRe, PAcIm, PAc);

            numericalDifferential numDiff
            (
                PAc,
                mesh,
                bubbleCloud.constProps().dict()
            );

            numDiff.firstOrderDiffWaveNumber
            (
                kSqrRe,
                kSqrIm,
                PCmpx,
                kSqrCloudInterp,
                dKidPr,
                dKrdPr,
                dKidPi,
                dKrdPi
            );

            // Build linear system with the jacobian
            volScalarField kSqrImPAcIm = kSqrIm + dKidPi*PAcIm - dKrdPi*PAcRe;
            volScalarField kSqrImPAcRe = kSqrIm + dKidPr*PAcRe + dKrdPr*PAcIm;
            fvScalarMatrix JacReEqn
            (
                fvm::laplacian(dPAcRe)
                + fvm::Sp(kSqrRe + dKrdPr*PAcRe - dKidPr*PAcIm, dPAcRe)
                - kSqrImPAcIm*dPAcIm
            );

            fvScalarMatrix JacImEqn
            (
                fvm::laplacian(dPAcIm)
                + fvm::Sp(kSqrRe + dKrdPi*PAcIm + dKidPi*PAcRe, dPAcIm)
                + kSqrImPAcRe*dPAcRe
            );


            // Prepare block matrix
            fvBlockMatrix<vector2> helmholtzEqnJac(blockJac);

            // Insert equations into block matrix
            helmholtzEqnJac.insertEquation(0, JacReEqn);
            helmholtzEqnJac.insertEquation(1, JacImEqn);

            // Add off-diagonal coupling terms
            helmholtzEqnJac.insertEquationCoupling(0, 1, -kSqrImPAcIm);
            helmholtzEqnJac.insertEquationCoupling(1, 0, kSqrImPAcRe);

            // Update source coupling: coupling terms eliminated from source
            helmholtzEqnJac.updateSourceCoupling();

            // Choose parallel or serial MUMPS and initialize
            autoPtr<blockMUMPSSolver> directSolver =
                    blockMUMPSSolver::New(helmholtzEqnJac, mesh);

            newtonRaphsonMethod nonLinSolver
            (
                loopControl,
                directSolver,
                helmholtzEqnJac,
                PAcRe,
                PAcIm,
                kSqrRe,
                kSqrIm,
                kSqrCloudInterp
            );

            nonLinSolver.solve(dPAcRe, dPAcIm);

            PAc = sqrt(sqr(PAcRe) + sqr(PAcIm));
        }


        //// Bjerknes force computation
        //
        volVectorField gradPAcRe = fvc::grad(PAcRe);
        volVectorField gradPAcIm = fvc::grad(PAcIm);

        // Workaround to apply atan2() on both internal and boundary fields.
        phiPAc.internalField() = Foam::atan2(PAcIm, PAcRe);
        phiPAc.boundaryField() =
            Foam::atan2(PAcIm.boundaryField(), PAcRe.boundaryField());

        absAndPhaseOfVectorField(G, psi, gradPAcRe, gradPAcIm);


        //// URANS
        // A combination of pimpleFoam (before consistency
        // update in foam-extend) and icoLagrangianFoam
        // PIMPLE loop with consistent formulation from pimpleDyMFoam
#       include "volContinuity.H"

        if (correctPhi && meshChanged)
        {
            // Fluxes will be corrected to absolute velocity
            // HJ, 6/Feb/2009
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

        if (meshChanged)
        {
#           include "CourantNo.H"
        }

        // --- PIMPLE loop
        while (pimple.loop())
        {
#           include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.H"
            }

            turbulence->correct();
        }

        runTime.write();

        scalar maxPAc = gMax(PAc);
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  Max(PAc) = " << maxPAc << " Pa" << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
