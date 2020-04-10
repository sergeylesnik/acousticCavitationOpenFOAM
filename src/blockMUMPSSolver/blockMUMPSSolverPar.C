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

#include "blockMUMPSSolverPar.H"
#include "addToRunTimeSelectionTable.H"
#include "IPstream.H"
#include "OPstream.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockMUMPSSolverPar, 0)
    addToRunTimeSelectionTable
    (
        blockMUMPSSolver,
        blockMUMPSSolverPar,
        runType
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMUMPSSolverPar::blockMUMPSSolverPar
(
    fvBlockMatrix<vector2> &matrix,
    const fvMesh &mesh
)
:
    blockMUMPSSolver(matrix, mesh),
    nNzLoc_(0),
    nRhsLoc_(0),
    nProcFaces_(0),
    localCellProcAddr_
    (
        IOobject
        (
        "cellProcAddressing",
        mesh_.facesInstance(),
        mesh_.meshSubDir,
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    ),
    lAddr_(matrix.lduAddr().lowerAddr()),
    uAddr_(matrix.lduAddr().upperAddr()),
    rhsGlobalIDs_(Pstream::nProcs()),
    dumpLinSysFromProcs_
    (
        MUMPSdict_.lookupOrDefault<Switch>("dumpLinSysFromProcs", 0)
    ),
    dumpProcIDs_
    (
        MUMPSdict_.lookupOrDefault<List<label> >
        (
            "dumpProcIDList", List<label>::zero
        )
    )
{
    irn_.resize(Pstream::nProcs());
    jcn_.resize(Pstream::nProcs());
    amv_.resize(Pstream::nProcs());
    rhs_.resize(Pstream::nProcs());
    myid_ = Pstream::myProcNo();

    exchangeProcPatches();

    // Matrix is non-symmetric for MUMPS since the diagonal
    // has coupling terms with different signs.
    // To optimize proc load only half of off-diagonal coeffs
    // known to the current proc is stored and sent to master.
    // Other half of coeffs is handled by the neigbouring proc.
    // Thus -> 2*nFaces entries in A for processor faces;
    // internal faces are handled completely by the current proc
    // -> 2 Entries * 2 Fluxes (upper&lower matrices) * nInternalFaces
    // diagonal is 2x2 -> 4*nCells entries in the matrix A for diagonal
    nNzLoc_ = 4*nCells_ + 4*nInternalFaces_ + 2*nProcFaces_;
    nRhsLoc_ = 2*nCells_;

    irn_[myid_].resize(nNzLoc_);
    jcn_[myid_].resize(nNzLoc_);
    amv_[myid_].resize(nNzLoc_);
    rhs_[myid_].resize(nRhsLoc_);
    rhsGlobalIDs_[myid_].resize(nRhsLoc_);

    assembleDiag();

    // Save the global ordering of the matrix rows for the corresp.
    // reordering of the rhs on the master later.
    // 4 entries per one cell in the diagonal block.
    // Every second entry in irn_ contains next row index.
    label rowI = 0;
    for (label i=0; i<(4*nCells_); i+=2)
    {
        rhsGlobalIDs_[myid_][rowI] = irn_[myid_][i];
        rowI++;
    }

    assembleOffDiag();
    assembleProcPatchOffDiag();

    mumps_.nz_loc = nNzLoc_;
    mumps_.irn_loc = irn_[myid_].data();
    mumps_.jcn_loc = jcn_[myid_].data();
    mumps_.a_loc = amv_[myid_].data();

    gatherProblemSize();

    // The input format of the distributed matrix
    mumps_.ICNTL(18) = 3;
    analyzeAndFactorizeMUMPS();
}


void Foam::blockMUMPSSolverPar::solveCore()
{
    if (debug)
    {
        Info<< "Gathering the right hand side on master before MUMPS call"
            << endl;
    }
    gatherRhs();

    solveMUMPS();

    if (debug)
    {
        Info<< "Redistributing solution to the slaves acc. to OF mapping"
            << endl;
    }
    scatterSol();
}


Foam::label Foam::blockMUMPSSolverPar::diagInd
(
    label i
)
{
    return 2 * localCellProcAddr_[i];
}

void Foam::blockMUMPSSolverPar::setOwnerNeighbInd
(
    label &ownerI,
    label &neighbourI,
    const label &faceI
)
{
    // Indices of the owners of the face
    ownerI = 2 * localCellProcAddr_[lAddr_[faceI]];

    // Indices of the neighbours of the owner cells
    // which share the same face
    neighbourI = 2 * localCellProcAddr_[uAddr_[faceI]];
}

void Foam::blockMUMPSSolverPar::exchangeProcPatches()
{
    label procPatchI = 0;
    nProcFaces_ = 0;

    // Cast procPolyPatches and put them in a list
    forAll(mesh_.boundary(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (isA<processorPolyPatch>(pp))
        {
            // Functions from processorPolyPatch needed.
            // Thus cast polyPatch to processorPolyPatch.
            const processorPolyPatch& procPP =
                  dynamicCast<const processorPolyPatch&>(pp);

            // Append ptr at the end of the PtrList.
            // Minimal size of List is 1.
            procPatchPtrs_.resize(procPatchI+1);
            procPatchPtrs_.set
            (
                procPatchI,
                new processorPolyPatch(procPP, procPP.boundaryMesh())
            );

            nProcFaces_ += procPP.size();
            procPatchIDs_.append(patchI);
            procPatchI++;
        }
    }

    if (debug)
    {
        Info<< "Exchanging cellProcAddressing between neighbouring processors"
            << endl;
    }

    neighbProcFaceOwners_.setSize(Pstream::nProcs());

    label nNeighbourFaces;

    // Only neigbour procs need to exchange. Set up own communication.
    forAll(procPatchPtrs_, patchI)
    {
        const int& neighbourProcNo = procPatchPtrs_[patchI].neighbProcNo();
        OPstream toNeighbour(Pstream::blocking, neighbourProcNo);
        toNeighbour << procPatchPtrs_[patchI].size();

        // Prepare array for sending in one send call
        List<label> sendBuf(procPatchPtrs_[patchI].size());

        label i=0;
        for
        (
            label faceI = procPatchPtrs_[patchI].start();
            faceI < procPatchPtrs_[patchI].start() +
                        procPatchPtrs_[patchI].size();
            faceI++
        )
        {
            sendBuf[i] = localCellProcAddr_[mesh_.faceOwner()[faceI]];
            i++;
        }

        toNeighbour << sendBuf;
    }

    forAll(procPatchPtrs_, patchI)
    {
        const int& neighbourProcNo = procPatchPtrs_[patchI].neighbProcNo();
        IPstream fromNeigbour(Pstream::blocking, neighbourProcNo);

        fromNeigbour >> nNeighbourFaces;

        // Assume that two neigbour procs share only one common patch
        neighbProcFaceOwners_[neighbourProcNo].setSize(nNeighbourFaces);

        fromNeigbour >> neighbProcFaceOwners_[neighbourProcNo];
    }
}

void Foam::blockMUMPSSolverPar::assembleProcPatchOffDiag()
{
    if (debug)
    {
        Info<< "Adding coeffs corresponding to the processor patches"
            << endl;
    }

    label ownerI;
    label neighbourI;

    forAll(procPatchPtrs_, patchI)
    {
        const int& neighbourProcNo = procPatchPtrs_[patchI].neighbProcNo();
        label patchGlobalI = procPatchIDs_[patchI];
        CoeffField<vector2>& patchCoeffField =
                matrix_.coupleUpper()[patchGlobalI];

        // Proc with the lower ID defines the upper matrix.
        // Proc with higher ID defines the lower matrix.
        // Interchange of owner and neigbour indices not needed
        // cause owner have automatically higher indices than the
        // neighbour from the other proc (compare to internal faces).
        for
        (
            label procFaceI = 0;
            procFaceI < procPatchPtrs_[patchI].size();
            procFaceI++
        )
        {
            // MUMPS; the indices of the C array start with 0
            // but the indices of the matrix start with 1 (acc. to MUMPS).
            // Consider only upper matrix since it's equal to the lower
            // and is not stored in BlockLduMatrix.
            // Coeffs of processor faces are multiplied by -1
            // -> premultiply by -1 to obtain true value

            // Order the local owner index to the global one
            ownerI = 2 * localCellProcAddr_
                    [
                        // Get the local owner index of the procFaceI
                        matrix_.lduAddr().patchAddr(patchGlobalI)[procFaceI]
                    ];

            neighbourI =
                    2 * neighbProcFaceOwners_[neighbourProcNo][procFaceI];

            addMatrixEntry
            (
                ownerI+1,
                neighbourI+1,
                -patchCoeffField.asLinear()[procFaceI][0]
            );
            addMatrixEntry
            (
                ownerI+2,
                neighbourI+2,
                -patchCoeffField.asLinear()[procFaceI][1]
            );
        }
    }
}

void Foam::blockMUMPSSolverPar::gatherProblemSize()
{

    label nNzGather = nNzLoc_;
    label nRhsGather = nRhsLoc_;
    Pstream::gather(nNzGather, sumOp<label>());
    Pstream::gather(nRhsGather, sumOp<label>());
    Pstream::gatherList(rhsGlobalIDs_);

    if (Pstream::master())
    {
        // Provide MUMPS nz and n variables on master only
        mumps_.nz = nNzGather;
        mumps_.n = nRhsGather;
        rhsMaster_.resize(nRhsGather);
    }
}

void Foam::blockMUMPSSolverPar::gatherRhs()
{
    Pstream::gatherList(rhs_);

    // rhsGlobalIDs_ hast to be gathered and rhsMaster_ resized before.
    // Order the RHS on Master acc. to the global row index.
    // Note that the first entry is 1 in MUMPS and 0 in OF List.
    if (Pstream::master())
    {
        forAll(rhs_, procI)
        {
            forAll(rhs_[procI], j)
            {
                rhsMaster_[rhsGlobalIDs_[procI][j]-1] = rhs_[procI][j];
            }
        }
        // Provide MUMPS the complete rhs on master
        mumps_.rhs = rhsMaster_.data();
    }
}

void Foam::blockMUMPSSolverPar::scatterSol()
{
    // rhs_ contains solution after MUMPS call.
    // Order the RHS on Slaves acc. to the OF mapping.
    // Note that the first entry is 1 in MUMPS and 0 in OF List.
    if (Pstream::master())
    {
        forAll(rhs_, procI)
        {
            forAll(rhs_[procI], j)
            {
                rhs_[procI][j] = rhsMaster_[rhsGlobalIDs_[procI][j]-1];
            }
        }
    }

    // scatterList() doesn't work. listCombineScatter() seems to do the job.
    Pstream::listCombineScatter(rhs_);
}

void Foam::blockMUMPSSolverPar::dumpCompleteRhs(fileName filePrefix)
{
    gatherProblemSize();
    gatherRhs();

    if (myid_ == 0)
    {
        OFstream osRhs
        (
            filePrefix + "_RHS_Time_" + mesh_.time().timeName()
        );

        osRhs << "procNo localI globRow rhs"<< endl;

        forAll(rhs_, procI)
        {
            forAll(rhs_[procI], i)
            {
                osRhs.precision(12);
                osRhs << procI << " " << i
                      << " " << rhsGlobalIDs_[procI][i]
                      << " " << rhs_[procI][i] << endl;
            }
        }
    }
}

void Foam::blockMUMPSSolverPar::dumpProcLinSys
(
    label procI,
    fileName filePrefix
)
{
    if (myid_ == procI)
    {
        OFstream osMatrix
        (
            filePrefix + "Matrix_Proc" + name(procI)
            + "_Time" + mesh_.time().timeName()
        );

        // Header
        osMatrix << "procNo localI globRow globCol coeff" << endl;

        forAll(irn_[procI], i)
        {
            osMatrix.precision(12);
            osMatrix << procI << " " << i
               << " " << irn_[procI][i]
               << " " << jcn_[procI][i]
               << " " << amv_[procI][i] << endl;
        }

        OFstream osRhs
        (
            filePrefix + "RHS_Proc" + name(procI)
            + "_Time" + mesh_.time().timeName()
        );

        // Header
        osRhs << "procNo localI globRow rhs"<< endl;

        forAll(rhs_[procI], i)
        {
            osRhs.precision(12);
            osRhs << procI << " " << i
                  << " " << rhsGlobalIDs_[procI][i]
                  << " " << rhs_[procI][i] << endl;
        }
    }
}

void Foam::blockMUMPSSolverPar::dumpAccToMUMPSDict()
{
    if (dumpCompleteLinSys_)
    {
        dumpMasterLinearSystem(dumpFilePrefix_);
    }

    if (dumpLinSysFromProcs_)
    {
        if (debug)
        {
            Info<< "Dumping proc linear system for the following procs: "
                << dumpProcIDs_ << endl;
        }

        forAll(dumpProcIDs_, procI)
        {
            dumpProcLinSys
            (
                dumpProcIDs_[procI],
                dumpFilePrefix_
            );
        }
    }
}
