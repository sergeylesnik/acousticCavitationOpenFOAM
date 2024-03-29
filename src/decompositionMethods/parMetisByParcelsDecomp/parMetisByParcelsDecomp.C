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

#include "parMetisByParcelsDecomp.H"
#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"
#include "labelIOField.H"
#include "syncTools.H"
#include "globalIndex.H"

#include <mpi.h>

extern "C"
{
#   include "parmetis.h"
}

#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parMetisByParcelsDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        parMetisByParcelsDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Does prevention of 0 cell domains and calls parmetis.
Foam::label Foam::parMetisByParcelsDecomp::decompose
(
    Field<label>& xadj,
    Field<label>& adjncy,
    const pointField& cellCentres,
    Field<label>& cellWeights,
    Field<label>& faceWeights,
    const labelList& options,
    labelList& finalDecomp
)
{
    // C style numbering
    label numFlag = 0;

    // Number of dimensions
    label nDims = 3;


    if (cellCentres.size() != xadj.size()-1)
    {
        FatalErrorIn("parMetisByParcelsDecomp::decompose(..)")
            << "cellCentres:" << cellCentres.size()
            << " xadj:" << xadj.size()
            << abort(FatalError);
    }


    // Get number of cells on all processors
    labelList nLocalCells(Pstream::nProcs());
    nLocalCells[Pstream::myProcNo()] = xadj.size()-1;
    Pstream::gatherList(nLocalCells);
    Pstream::scatterList(nLocalCells);

    // Get cell offsets.
    labelList cellOffsets(Pstream::nProcs()+1);
    label nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }
    cellOffsets[Pstream::nProcs()] = nGlobalCells;

    // Convert pointField into the data type parMetis expects (float or double)
    Field<real_t> xyz(3*cellCentres.size());
    label compI = 0;
    forAll(cellCentres, cellI)
    {
        const point& cc = cellCentres[cellI];
        xyz[compI++] = float(cc.x());
        xyz[compI++] = float(cc.y());
        xyz[compI++] = float(cc.z());
    }

    // Make sure every domain has at least one cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (Metis falls over with zero sized domains)
    // Trickle cells from processors that have them up to those that
    // don't.


    // Number of cells to send to the next processor
    // (is same as number of cells next processor has to receive)
    labelList nSendCells(Pstream::nProcs(), 0);

    for (label procI = nLocalCells.size()-1; procI >=1; procI--)
    {
        if (nLocalCells[procI]-nSendCells[procI] < 1)
        {
            nSendCells[procI-1] = nSendCells[procI]-nLocalCells[procI]+1;
        }
    }

    // First receive (so increasing the sizes of all arrays)

    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        // Receive cells from previous processor
        IPstream fromPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        Field<label> prevXadj(fromPrevProc);
        Field<label> prevAdjncy(fromPrevProc);
        Field<real_t> prevXyz(fromPrevProc);
        Field<label> prevCellWeights(fromPrevProc);
        Field<label> prevFaceWeights(fromPrevProc);

        if (prevXadj.size() != nSendCells[Pstream::myProcNo()-1])
        {
            FatalErrorIn("parMetisByParcelsDecomp::decompose(..)")
                << "Expected from processor " << Pstream::myProcNo()-1
                << " connectivity for " << nSendCells[Pstream::myProcNo()-1]
                << " nCells but only received " << prevXadj.size()
                << abort(FatalError);
        }

        // Insert adjncy
        prepend(prevAdjncy, adjncy);
        // Adapt offsets and prepend xadj
        xadj += prevAdjncy.size();
        prepend(prevXadj, xadj);
        // Coords
        prepend(prevXyz, xyz);
        // Weights
        prepend(prevCellWeights, cellWeights);
        prepend(prevFaceWeights, faceWeights);
    }


    // Send to my next processor

    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        // Send cells to next processor
        OPstream toNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        label nCells = nSendCells[Pstream::myProcNo()];
        label startCell = xadj.size()-1 - nCells;
        label startFace = xadj[startCell];
        label nFaces = adjncy.size()-startFace;

        // Send for all cell data: last nCells elements
        // Send for all face data: last nFaces elements
        toNextProc
            << Field<label>::subField(xadj, nCells, startCell)-startFace
            << Field<label>::subField(adjncy, nFaces, startFace)
            << SubField<real_t>(xyz, nDims*nCells, nDims*startCell)
            <<
            (
                cellWeights.size()
              ? static_cast<const Field<label>&>
                (
                    Field<label>::subField(cellWeights, nCells, startCell)
                )
              : Field<label>(0)
            )
            <<
            (
                faceWeights.size()
              ? static_cast<const Field<label>&>
                (
                    Field<label>::subField(faceWeights, nFaces, startFace)
                )
              : Field<label>(0)
            );

        // Remove data that has been sent
        if (faceWeights.size())
        {
            faceWeights.setSize(faceWeights.size()-nFaces);
        }
        if (cellWeights.size())
        {
            cellWeights.setSize(cellWeights.size()-nCells);
        }
        xyz.setSize(xyz.size()-nDims*nCells);
        adjncy.setSize(adjncy.size()-nFaces);
        xadj.setSize(xadj.size() - nCells);
    }



    // Adapt number of cells
    forAll(nSendCells, procI)
    {
        // Sent cells
        nLocalCells[procI] -= nSendCells[procI];

        if (procI >= 1)
        {
            // Received cells
            nLocalCells[procI] += nSendCells[procI-1];
        }
    }
    // Adapt cellOffsets
    nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }


    if (nLocalCells[Pstream::myProcNo()] != (xadj.size()-1))
    {
        FatalErrorIn("parMetisByParcelsDecomp::decompose(..)")
            << "Have connectivity for " << xadj.size()-1
            << " cells but nLocalCells:" << nLocalCells[Pstream::myProcNo()]
            << abort(FatalError);
    }

    // Weight info
    label wgtFlag = 0;
    label* vwgtPtr = nullptr;
    label* adjwgtPtr = nullptr;

    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
        wgtFlag += 2;       // Weights on vertices
    }
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
        wgtFlag += 1;       // Weights on edges
    }


    // Number of weights or balance constraints
    label nCon = 1;
    // Per processor, per constraint the weight
    Field<real_t> tpwgts(nCon*nProcessors_, 1./nProcessors_);
    // Imbalance tolerance
    Field<real_t> ubvec(nCon, 1.02);
    if (nProcessors_ == 1)
    {
        // If only one processor there is no imbalance.
        ubvec[0] = 1;
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    // output: cell -> processor addressing
    finalDecomp.setSize(nLocalCells[Pstream::myProcNo()]);

    // output: number of cut edges
    label edgeCut = 0;

    // Number of parts
    label nProcs = nProcessors_;
    
    // Ratio of inter-processor communication time compared to data
    // redistribution time. 1000.0 is advised by documentation.
    float itr = 1000.0;

    // This algorithm tries to preserve the most of the current local part of
    // the decomposition minimizing migration cost and MPI transfer data amount

    Info<< "Rebalancing mesh adaptively using current decomposition" << endl;
    ParMETIS_V3_AdaptiveRepart
    (
        cellOffsets.begin(),    // vtxDist
        xadj.begin(),
        adjncy.begin(),
        vwgtPtr,                // vertexweights
        NULL,
        adjwgtPtr,              // edgeweights
        &wgtFlag,
        &numFlag,
        &nCon,
        &nProcs,                // nParts
        tpwgts.begin(),
        ubvec.begin(),
        &itr,
        const_cast<List<label>&>(options).begin(),
        &edgeCut,
        finalDecomp.begin(),
        &comm
    );

    // Info<< "Rebalancing mesh without using current decomposition" << endl;
    // ParMETIS_V3_PartGeomKway
    // (
    //     cellOffsets.begin(),    // vtxDist
    //     xadj.begin(),
    //     adjncy.begin(),
    //     vwgtPtr,                // vertexweights
    //     adjwgtPtr,              // edgeweights
    //     &wgtFlag,
    //     &numFlag,
    //     &nDims,
    //     xyz.begin(),
    //     &nCon,
    //     &nProcs,                // nParts
    //     tpwgts.begin(),
    //     ubvec.begin(),
    //     const_cast<List<label>&>(options).begin(),
    //     &edgeCut,
    //     finalDecomp.begin(),
    //     &comm
    // );


    // If we sent cells across make sure we undo it
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Receive back from next processor if I sent something
    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        IPstream fromNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        labelList nextFinalDecomp(fromNextProc);

        if (nextFinalDecomp.size() != nSendCells[Pstream::myProcNo()])
        {
            FatalErrorIn("parMetisByParcelsDecomp::decompose(..)")
                << "Expected from processor " << Pstream::myProcNo()+1
                << " decomposition for " << nSendCells[Pstream::myProcNo()]
                << " nCells but only received " << nextFinalDecomp.size()
                << abort(FatalError);
        }

        append(nextFinalDecomp, finalDecomp);
    }

    // Send back to previous processor.
    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        OPstream toPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        label nToPrevious = nSendCells[Pstream::myProcNo()-1];

        toPrevProc <<
            SubList<label>
            (
                finalDecomp,
                nToPrevious,
                finalDecomp.size()-nToPrevious
            );

        // Remove locally what has been sent
        finalDecomp.setSize(finalDecomp.size()-nToPrevious);
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parMetisByParcelsDecomp::parMetisByParcelsDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::parMetisByParcelsDecomp::decompose
(
    const pointField& cc,
    const scalarField& cWeights
)
{
    if (cc.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "parMetisByParcelsDecomp::decompose"
            "(const pointField&, const scalarField&)"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << cc.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // For running sequential ...
    if (Pstream::nProcs() <= 1)
    {
        return metisDecomp(decompositionDict_, mesh_).decompose
        (
            cc,
            cWeights
        );
    }


    // Connections
    Field<label> adjncy;
    // Offsets into adjncy
    Field<label> xadj;
    calcDistributedCSR
    (
        mesh_,
        adjncy,
        xadj
    );


    // decomposition options. 0 = use defaults
    labelList options(3, label(0));
    //options[0] = 1;     // don't use defaults but use values below
    //options[1] = -1;    // full debug info
    //options[2] = 15;    // random number seed

    // cell weights (so on the vertices of the dual)
    Field<label> cellWeights;

    // face weights (so on the edges of the dual)
    Field<label> faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != mesh_.nCells())
        {
            FatalErrorIn
            (
                "parMetisByParcelsDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << mesh_.nCells()
                << exit(FatalError);
        }

        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("parMetisByParcelsCoeffs"))
    {
        const dictionary& parMetisByParcelsCoeffs =
            decompositionDict_.subDict("parMetisByParcelsCoeffs");
        word weightsFile;

        if (parMetisByParcelsCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            if (mesh_.found(weightsFile))
            {
                Info<< "parMetisByParcelsDecomp : Using an IOregistry object "
                    << "with cell-based weights: " << weightsFile << endl;
                labelField cellIOWeights
                (
                    mesh_.lookupObject<labelField>(weightsFile)
                );

                cellWeights.transfer(cellIOWeights);
            }
            else // no data in registry; read from disk
            {
                Info<< "parMetisByParcelsDecomp : Using cell-based weights read "
                    << "from " << weightsFile << endl;

                labelIOField cellIOWeights
                (
                    IOobject
                    (
                        weightsFile,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                );

                cellWeights.transfer(cellIOWeights);
            }

            if (cellWeights.size() != mesh_.nCells())
            {
                FatalErrorIn
                (
                    "parMetisByParcelsDecomp::decompose"
                    "(const pointField&, const scalarField&)"
                )   << "Number of cell weights " << cellWeights.size()
                    << " does not equal number of cells " << mesh_.nCells()
                    << exit(FatalError);
            }
        }

        if (parMetisByParcelsCoeffs.readIfPresent("faceWeightsFile", weightsFile))
        {
            Info<< "parMetisByParcelsDecomp : Using face-based weights read from "
                << weightsFile << endl;

            labelIOField weights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );

            if (weights.size() != mesh_.nFaces())
            {
                FatalErrorIn
                (
                    "parMetisByParcelsDecomp::decompose"
                    "(const pointField&, const scalarField&)"
                )   << "Number of face weights " << weights.size()
                    << " does not equal number of internal and boundary faces "
                    << mesh_.nFaces()
                    << exit(FatalError);
            }

            faceWeights.setSize(adjncy.size());

            // Assume symmetric weights. Keep same ordering as adjncy.
            labelList nFacesPerCell(mesh_.nCells(), 0);

            // Handle internal faces
            for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
            {
                label w = weights[faceI];

                label own = mesh_.faceOwner()[faceI];
                label nei = mesh_.faceNeighbour()[faceI];

                faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
                faceWeights[xadj[nei] + nFacesPerCell[nei]++] = w;
            }
            // Coupled boundary faces
            const polyBoundaryMesh& patches = mesh_.boundaryMesh();

            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (pp.coupled())
                {
                    label faceI = pp.start();

                    forAll(pp, i)
                    {
                        label w = weights[faceI];
                        label own = mesh_.faceOwner()[faceI];
                        faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
                        faceI++;
                    }
                }
            }
        }

        if (parMetisByParcelsCoeffs.readIfPresent("options", options))
        {
            Info<< "Using Metis options     " << options
                << nl << endl;

            if (options.size() != 3)
            {
                FatalErrorIn
                (
                    "parMetisByParcelsDecomp::decompose"
                    "(const pointField&, const scalarField&)"
                )   << "Number of options " << options.size()
                    << " should be three." << exit(FatalError);
            }
        }
    }


    // Do actual decomposition
    labelList finalDecomp;
    decompose
    (
        xadj,
        adjncy,
        cc,
        cellWeights,
        faceWeights,
        options,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


Foam::labelList Foam::parMetisByParcelsDecomp::decompose
(
    const labelList& cellToRegion,
    const pointField& regionPoints,
    const scalarField& regionWeights
)
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (cellToRegion.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "parMetisByParcelsDecomp::decompose(const labelList&, const pointField&)"
        )   << "Size of cell-to-coarse map " << cellToRegion.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }


    // Global region numbering engine
    globalIndex globalRegions(regionPoints.size());


    // Get renumbered owner region on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalNeighbour(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                label ownRegion = cellToRegion[faceOwner[faceI]];
                globalNeighbour[bFaceI++] = globalRegions.toGlobal(ownRegion);
                faceI++;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh_, globalNeighbour, false);


    // Get globalCellCells on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList globalRegionRegions;
    {
        List<dynamicLabelList > dynRegionRegions(regionPoints.size());

        // Internal faces first
        forAll(faceNeighbour, faceI)
        {
            label ownRegion = cellToRegion[faceOwner[faceI]];
            label neiRegion = cellToRegion[faceNeighbour[faceI]];

            if (ownRegion != neiRegion)
            {
                label globalOwn = globalRegions.toGlobal(ownRegion);
                label globalNei = globalRegions.toGlobal(neiRegion);

                if (findIndex(dynRegionRegions[ownRegion], globalNei) == -1)
                {
                    dynRegionRegions[ownRegion].append(globalNei);
                }
                if (findIndex(dynRegionRegions[neiRegion], globalOwn) == -1)
                {
                    dynRegionRegions[neiRegion].append(globalOwn);
                }
            }
        }

        // Coupled boundary faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                label faceI = pp.start();
                label bFaceI = pp.start() - mesh_.nInternalFaces();

                forAll(pp, i)
                {
                    label ownRegion = cellToRegion[faceOwner[faceI]];
                    label globalNei = globalNeighbour[bFaceI++];
                    faceI++;

                    if
                    (
                        findIndex(dynRegionRegions[ownRegion], globalNei)
                     == -1
                    )
                    {
                        dynRegionRegions[ownRegion].append(globalNei);
                    }
                }
            }
        }

        globalRegionRegions.setSize(dynRegionRegions.size());
        forAll(dynRegionRegions, i)
        {
            globalRegionRegions[i].transfer(dynRegionRegions[i]);
        }
    }

    labelList regionDecomp
    (
        decompose
        (
            globalRegionRegions,
            regionPoints,
            regionWeights
        )
    );

    // Rework back into decomposition for original mesh_
    labelList cellDistribution(cellToRegion.size());

    forAll(cellDistribution, cellI)
    {
        cellDistribution[cellI] = regionDecomp[cellToRegion[cellI]];
    }

    fixCyclics(mesh_, cellDistribution);

    return cellDistribution;
}


Foam::labelList Foam::parMetisByParcelsDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "parMetisByParcelsDecomp::decompose(const labelListList&"
            ", const pointField&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ") or weights (" << cWeights.size() << ")." << exit(FatalError);
    }

    // For running sequential ...
    if (Pstream::nProcs() <= 1)
    {
        return metisDecomp(decompositionDict_, mesh_).decompose
        (
            globalCellCells,
            cellCentres,
            cWeights
        );
    }


    // Make Metis Distributed CSR (Compressed Storage Format) storage

    // Connections
    Field<label> adjncy;

    // Offsets into adjncy
    Field<label> xadj;

    calcCSR(globalCellCells, adjncy, xadj);

    // decomposition options. 0 = use defaults
    labelList options(3, label(0));
    //options[0] = 1;     // don't use defaults but use values below
    //options[1] = -1;    // full debug info
    //options[2] = 15;    // random number seed

    // cell weights (so on the vertices of the dual)
    Field<label> cellWeights;

    // face weights (so on the edges of the dual)
    Field<label> faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "parMetisByParcelsDecomp::decompose(const labelListList&"
                ", const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != globalCellCells.size())
        {
            FatalErrorIn
            (
                "parMetisByParcelsDecomp::decompose(const labelListList&"
                ", const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << globalCellCells.size()
                << exit(FatalError);
        }

        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }


    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("parMetisByParcelsCoeffs"))
    {
        const dictionary& parMetisByParcelsCoeffs =
            decompositionDict_.subDict("parMetisByParcelsCoeffs");

        if (parMetisByParcelsCoeffs.readIfPresent("options", options))
        {
            Info<< "Using Metis options     " << options
                << nl << endl;

            if (options.size() != 3)
            {
                FatalErrorIn
                (
                    "parMetisByParcelsDecomp::decompose(const labelListList&"
                    ", const pointField&, const scalarField&)"
                )   << "Number of options " << options.size()
                    << " should be three." << exit(FatalError);
            }
        }
    }


    // Do actual decomposition
    labelList finalDecomp;
    decompose
    (
        xadj,
        adjncy,
        cellCentres,
        cellWeights,
        faceWeights,
        options,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


// ************************************************************************* //
