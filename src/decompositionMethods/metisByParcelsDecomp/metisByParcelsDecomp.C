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

#include "metisByParcelsDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"

#include "volFields.H"
#include "fvCFD.H"


extern "C"
{
#define OMPI_SKIP_MPICXX
#   include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisByParcelsDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisByParcelsDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisByParcelsDecomp::decompose
(
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& cWeights,
    labelList& finalDecomp
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("k-way");

    label numCells = xadj.size()-1;

    // decomposition options. 0 = use defaults
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // processor weights initialised with no size, only used if specified in
    // a file - use the data type that metis expects here
    Field<real_t> processorWeights;

    // cell weights (so on the vertices of the dual)
    labelList cellWeights;

    // face weights (so on the edges of the dual)
    labelList faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "metisByParcelsDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorIn
            (
                "metisByParcelsDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
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
    if (decompositionDict_.found("metisByParcelsCoeffs"))
    {
        const dictionary& metisByParcelsCoeffs =
            decompositionDict_.subDict("metisByParcelsCoeffs");
        word weightsFile;

        if (metisByParcelsCoeffs.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorIn("metisByParcelsDecomp::decompose()")
                    << "Method " << method
                    << " in metisByParcelsCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisByParcelsDecomp : Using Metis method     " << method
                << nl << endl;
        }

        List<int> mOptions;
        if (metisByParcelsCoeffs.readIfPresent("options", mOptions))
        {
            if (mOptions.size() != METIS_NOPTIONS)
            {
                FatalErrorIn("metisByParcelsDecomp::decompose()")
                    << "Number of options in metisByParcelsCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            forAll(mOptions, i)
            {
                options[i] = mOptions[i];
            }

            Info<< "metisByParcelsDecomp : Using Metis options     " << mOptions
                << nl << endl;
        }

        if (metisByParcelsCoeffs.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != nProcessors_)
            {
                FatalErrorIn("metisByParcelsDecomp::decompose(const pointField&)")
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nProcessors_
                    << exit(FatalError);
            }
        }

        if (metisByParcelsCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "metisByParcelsDecomp : Using cell-based weights." << endl;

            // Create a mesh object in order to work with volScalarField
            Foam::fvMesh mesh
            (
                Foam::IOobject
                (
                    Foam::fvMesh::defaultRegion,
                    mesh_.time().timeName(),
                    mesh_.time(),
                    Foam::IOobject::MUST_READ
                )
            );

            volScalarField volScalarCellWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                1
            );

            if (volScalarCellWeights.headerOk())
            {
                Info<< "Cell weights file found and will be used for "
                    << "processor balancing: " << weightsFile
                    << endl;
            }
            else
            {
                Info<< "No cell weights file found. Same weights will be used "
                    << "for all cells." << endl;
            }
            
            // Convert scalarField to labelField
            labelField labelCellWeights(mesh_.nCells());
            forAll(labelCellWeights, cellI)
            {
                labelCellWeights[cellI] = volScalarCellWeights.field()[cellI];
            }

            cellWeights.transfer(labelCellWeights);

            if (cellWeights.size() != xadj.size()-1)
            {
                FatalErrorIn("metisByParcelsDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " does not equal number of cells " << xadj.size()-1
                    << exit(FatalError);
            }
        }
    }

    label nProcs = nProcessors_;

    // output: cell -> processor addressing
    finalDecomp.setSize(numCells);

    // output: number of cut edges
    label edgeCut = 0;

    // Vertex weight info
    label* vwgtPtr = nullptr;
    label* adjwgtPtr = nullptr;

    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
    }
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
    }

    label one = 1;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<label>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<label>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            nullptr,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            nullptr,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<label>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<label>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            nullptr,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            nullptr,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisByParcelsDecomp::metisByParcelsDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisByParcelsDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "metisByParcelsDecomp::decompose(const pointField&,const scalarField&)"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    labelList adjncy;
    labelList xadj;
    calcCSR
    (
        mesh_,
        adjncy,
        xadj
    );

    // Decompose using default weights
    labelList finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


Foam::labelList Foam::metisByParcelsDecomp::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    if (fineToCoarse.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "metisByParcelsDecomp::decompose"
            "(const labelList&, const pointField&, const scalarField&)"
        )   << "Size of cell-to-coarse map " << fineToCoarse.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    labelList adjncy;
    labelList xadj;
    {
        // Get cellCells on coarse mesh.
        labelListList cellCells;
        calcCellCells
        (
            mesh_,
            fineToCoarse,
            coarsePoints.size(),
            cellCells
        );

        calcCSR(cellCells, adjncy, xadj);
    }

    // Decompose using default weights
    labelList finalDecomp;
    decompose(adjncy, xadj, coarseWeights, finalDecomp);


    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[fineToCoarse[i]];
    }

    fixCyclics(mesh_, fineDistribution);

    return fineDistribution;
}


Foam::labelList Foam::metisByParcelsDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights
)
{
    if (cc.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "metisByParcelsDecomp::decompose"
            "(const pointField&, const labelListList&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cc.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    labelList adjncy;
    labelList xadj;
    calcCSR(globalCellCells, adjncy, xadj);

    // Decompose using default weights
    labelList finalDecomp;
    decompose(adjncy, xadj, cWeights, finalDecomp);

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
