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
    blockMUMPSSolverPar

Description
    An implementation of the abstract blockMUMPSSolver for the solution of
    linear systems with the direct solver MUMPS in parallel.

Author
    Sergey Lesnik, ITM Clausthal

SourceFiles
    blockMUMPSSolverPar.C

\*---------------------------------------------------------------------------*/

#ifndef blockMUMPSSolverPar_H
#define blockMUMPSSolverPar_H

#include "blockMUMPSSolver.H"
#include "typeInfo.H"
#include "dynamicLabelList.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class blockMUMPSSolverPar Declaration
\*---------------------------------------------------------------------------*/

class blockMUMPSSolverPar
:
    public blockMUMPSSolver
{
    // Private data

        //- Number of non-zeros in the matrix corresp. to the proc
        label nNzLoc_;

        //- Number of RHS elements on the proc
        label nRhsLoc_;

        //- Number of proc-to-proc faces on the proc
        label nProcFaces_;

        //- Addressing of local cells on the proc in the global scope
        labelIOList localCellProcAddr_;

        //- Internal faces on the processor: lower and upper addressing
        const labelList& lAddr_;
        const labelList& uAddr_;

        //- List of patches neigbouring to other processors
        PtrList<processorPolyPatch> procPatchPtrs_;

        //- IDs of the patches for retrieving of matrix coeffs
        dynamicLabelList procPatchIDs_;

        //- IDs of the face owners from the neigbour proc
        labelListList neighbProcFaceOwners_;

        //- Global IDs of rows on the proc needed for RHS reordering on master
        labelListList rhsGlobalIDs_;

        //- Right hand side of the whole linear system; used on master only.
        //  Solution is stored after the MUMPS solver call.
        scalarList rhsMaster_;


        //
        const bool dumpLinSysFromProcs_;
        labelList dumpProcIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blockMUMPSSolverPar(const blockMUMPSSolverPar&);

        //- Disallow default bitwise assignment
        void operator=(const blockMUMPSSolverPar&);


public:

    //- Runtime type information
    TypeName("blockMUMPSSolverPar")

    // Constructors

        //- Construct from matrix components and mesh
        blockMUMPSSolverPar
        (
            fvBlockMatrix<vector2>& matrix,
            const fvMesh& mesh
        );

    // Destructor

        virtual ~blockMUMPSSolverPar()
        {}

    // Member Functions

        //- Solve with different RHS from provided matrix.
        //  Useful when LHS doesn't change.
        void solveCore();

        //- Get global index for MUMPS of the diagonal entries
        label diagInd(label i);

        //- Set global indices for MUMPS of the owner and neighbour for a face
        void setOwnerNeighbInd
        (
            label &ownerI,
            label &neighbourI,
            const label &faceI
        );

        //- Exchange proc patches between neighbour procs
        void exchangeProcPatches();

        //- Assemble off-diagonal terms on proc patches
        void assembleProcPatchOffDiag();

        //- Gather matrix orders, RHS global IDs and sizes from slaves
        void gatherProblemSize();

        //- Gather RHS on master and order it acc. to glabl IDs
        void gatherRhs();

        //- Redistribute solution from master to slaves
        void scatterSol();

        //- Gather and dump current complete RHS
        void dumpCompleteRhs(fileName filePrefix);

        //- Dump partial linear system (matrix and RHS) from a specified proc
        void dumpProcLinSys(label procI, fileName filePrefix);

        //- Dump according to MUMPS dictionary from case
        void dumpAccToMUMPSDict();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // blockMUMPSSolverPar_H

// ************************************************************************* //
