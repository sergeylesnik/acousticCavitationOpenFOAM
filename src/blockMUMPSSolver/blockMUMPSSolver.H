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

Class
    blockMUMPSSolver

Description
    An interface class for the usage of the direct linear solver MUMPS.

Author
    Sergey Lesnik, ITM Clausthal

SourceFiles
    blockMUMPSSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blockMUMPSSolver_H
#define blockMUMPSSolver_H

// fvc.H is not required by this class but it has to be included because it's
// erroneously absent in fvBlockMatrix.H
#include "fvc.H"

#include "parRun.H"
#include "fvBlockMatrix.H"
#include "vector2.H"
#include "blockLduSolvers.H"
#include "labelIOList.H"
#include "fvMesh.H"
#include "foamTime.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "Pstream.H"

#include "mpi.h"

#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

// Macro substituting indices in order to match MUMPS documentation
#define ICNTL(I) icntl[(I)-1]

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blockMUMPSSolver Declaration
\*---------------------------------------------------------------------------*/

class blockMUMPSSolver

{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        blockMUMPSSolver(const blockMUMPSSolver&);

        //- Disallow default bitwise assignment
        void operator=(const blockMUMPSSolver&);


protected:

        //- Matrix&Mesh
        fvBlockMatrix<vector2>& matrix_;
        const fvMesh& mesh_;

        //- Number of cells and internal faces per proc
        const label nCells_;
        const label nInternalFaces_;

        //- Arrays of arrays which hold the matrix data ordered by proc Nr.
        labelListList irn_;  // array of row indices
        labelListList jcn_;  // array of column indices
        scalarListList amv_;  // array of values of the matrix A
        scalarListList rhs_;  // array of values of the RHS

        //- Processor ID
        label myid_;

        //- Index of current matrix entry
        label matrixI_;

        // MUMPS instance
        DMUMPS_STRUC_C mumps_;

        // Normalized residuals
        vector2 initialNormResidual_;
        vector2 finalNormResidual_;

        // MUMPS dictionary and its entries from case dir
        const IOdictionary& MUMPSdict_;
        const bool dumpCompleteLinSys_;
        const word dumpFilePrefix_;
        const bool printResiduals_;


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockMUMPSSolver,
            runType,
            (
                fvBlockMatrix<vector2> &matrix,
                const fvMesh& mesh
            ),
            (
                matrix,
                mesh
            )
        );

    TypeName("blockMUMPSSolver");


    // Constructors

        //- Construct from matrix components and mesh
        blockMUMPSSolver
        (
            fvBlockMatrix<vector2>& matrix,
            const fvMesh& mesh
        );


    // Selectors

        //- OF factory method to choose appropriate derived class.
        //  In this case: serial or parallel version of MUMPS.
        static autoPtr<blockMUMPSSolver> New
        (
            fvBlockMatrix<vector2> &matrix,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~blockMUMPSSolver()
        {
            // Terminate MUMPS instance
            mumps_.job = JOB_END;
            dmumps_c(&mumps_);
        }


    // Member Functions

        //- Solve routine
        void solve();

        //- Solve with different RHS from the provided matrix
        //  Useful when LHS doesn't change.
        void solveWithRhs(const fvBlockMatrix<vector2>& matrix);

        virtual void solveCore() = 0;

        //- Compute and store residual
        //  Computation procedure is OF alike, but ignores boundary and
        //  interface faces.
        vector2 computeResidual
        (
            const Field<vector2>& x,
            const Field<vector2> &b,
            const fvBlockMatrix<vector2>& matrix
        );

        //- Return initial residual
        const vector2& initialResidual() const
        {
            return initialNormResidual_;
        }

        //- Print normalized initial and final residuals
        void printResiduals();

        //- Analyze and factorize provided matrix by MUMPS.
        //  RHS doesn't need to be defined.
        void analyzeAndFactorizeMUMPS();

        //- Set up and execute MUMPS
        void solveMUMPS();

        //- Assemble diagonal
        void assembleDiag();

        //- Assemble upper and lower matrices
        void assembleOffDiag();

        //- Get global index for MUMPS of the diagonal entries
        virtual label diagInd(label i) = 0;

        //- Set global indices for MUMPS of the owner and neighbour for a face
        virtual void setOwnerNeighbInd
        (
            label &ownerI,
            label &neighbourI,
            const label &faceI
        ) = 0;

        //- Add row, column number and corresp. coefficient to class arrays
        void addMatrixEntry(label rowI, label colI, scalar aCoeff);

        //- Get RHS from OF linear system class (source())
        void getRhs(const fvBlockMatrix<vector2>& matrix);

        void getRhs();

        //- Pass solution to OF linear system class (psi()). rhs_ has to
        //  contain solution which is the case after MUMPS execution
        void passSol();

        //- Dump current linear system set on master to hard-drive
        void dumpMasterLinearSystem(fileName filePrefix);

        //- Dump current complete RHS
        virtual void dumpCompleteRhs(fileName filePrefix) = 0;

        //- Dump according to MUMPS dictionary from case
        virtual void dumpAccToMUMPSDict() = 0;

    // Access

        inline const fvMesh& mesh() const
        {
            return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // blockMUMPSSolver_H

// ************************************************************************* //
