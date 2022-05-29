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

#include "loadBalanceChecker.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::loadBalanceChecker::loadBalanceChecker
(
    const fvMesh &mesh,
    const AcCavitationCloud<acCavitationParcel> &cloud
)
:
    mesh_(mesh),
    cloud_(cloud),
    dynamicMeshDict_
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh_.time().caseConstant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    loadBalanceCoeffs_(dynamicMeshDict_.subDict("myLoadBalanceFvMeshCoeffs")),
    imbalanceTrigger_
    (
        readScalar(loadBalanceCoeffs_.lookup("imbalanceTrigger"))
    ),
    parMetisByParcelsCoeffs_
    (
        loadBalanceCoeffs_.subDict("parMetisByParcelsCoeffs")
    ),
    weightsFileName_(parMetisByParcelsCoeffs_.lookup("cellWeightsFile")),
    cellIOWeights_
    (
        IOobject
        (
            weightsFileName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_.nCells()
    ),
    injectedParcelsWeightAmplifier_
    (
        readScalar
        (
            loadBalanceCoeffs_.lookup("injectedParcelsWeightAmplifier")
        )
    )
{
    const word& decomposeMethod = loadBalanceCoeffs_.lookup("method");

    if (decomposeMethod != "parMetisByParcels")
    {
        FatalErrorIn
        (
            "loadBalanceChecker::loadBalanceChecker()"
        )   << "Load balancing works only with the decomposition method "
            << "parMetisByParcels but " << decomposeMethod
            << " is specified in " << dynamicMeshDict_.filePath()
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loadBalanceChecker::isImba()
{
    if (Pstream::parRun())
    {
        labelList parcelsPerCell = cloud_.nParcelsPerCell();

        // Mesh size of the proc might have changed after load balancing
        cellIOWeights_.resize(mesh_.nCells());

        forAll(parcelsPerCell, cellI)
        {
            // Make sure there are no zeros in the list. Metis doesn't like it.
            cellIOWeights_[cellI] = max(1, parcelsPerCell[cellI]);
        }

        // Calculate local and global load
        scalar localLoad = sum(cellIOWeights_);
        scalar minWeights = returnReduce(localLoad, minOp<scalar>());
        scalar maxWeights = returnReduce(localLoad, maxOp<scalar>());
        scalar injectedParcelsLoad =
            injectedParcelsWeightAmplifier_
          * (
                cloud_.parcelsAddedConti()
              + cloud_.injection().parcelsAddedTotal()
            );
        localLoad += injectedParcelsLoad;

        scalar globalLoad = localLoad;

        reduce(globalLoad, sumOp<scalar>());

        globalLoad /= Pstream::nProcs();

        // Calculate imbalance as min of localLoad/globalLoad
        scalar imbalance = mag(1 - localLoad / globalLoad);

        reduce(imbalance, maxOp<scalar>());

        Info<< "Global load imbalance: " << imbalance
            << "  Average load: " << globalLoad
            << "  Number of parcels on a proc Min: " << minWeights
            << ", Max: " << maxWeights << endl;

        if (imbalance > imbalanceTrigger_)
        {
            return true;
        }
    }

    return false;
}
