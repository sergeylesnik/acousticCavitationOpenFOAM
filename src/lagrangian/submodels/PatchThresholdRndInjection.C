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

#include "PatchThresholdRndInjection.H"
#include "DataEntry.H"
#include "pdf.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchThresholdRndInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    notImplemented("PatchThresholdRndInjection<CloudType>::parcelsToInject");
    // time0 and time1 are relative to SOI
    return 0;
}


template<class CloudType>
Foam::scalar Foam::PatchThresholdRndInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    notImplemented("PatchThresholdRndInjection<CloudType>::volumeToInject");
    return 0.0;
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::prepareForNextTimeStep
(
    const scalar time,
    label &newParcels,
    scalar &newVolume
)
{
    // Calculate the number of parcels to be such that given void fraction
    // is not exceeded.

    newVolume = 0.0;
    newParcels = 0;
    currLocalCellI_ = 0;
    parcelDiameters_.clearStorage();

    if (!isInjectionTime())
    {
        return;
    }

    forAll(cellOwners_, cellI)
    {
        parcelsInjPerCellList_[cellI] = 0;

        label cellOwner = cellOwners_[cellI];
        scalar nPCurrCell = nParticlesPerParcel_[cellI];
        referenceFieldList_[cellI] = referenceField_[cellOwner];
        scalar voidFracI = referenceFieldList_[cellI];
        scalar cellVolume = this->owner().mesh().cellVolumes()[cellOwner];

        // Volume to inject is restricted by the threshold, current cell
        // void fraction. Also substract parcel volume of the PDF's median
        // diameter to fulfill threshold on average.
        scalar volumeInjCurrCell =
            (threshold_ - voidFracI) * cellVolume
            - nPCurrCell*medianParticleVolume_;

        // Restrict to positive volume
        volumeInjCurrCell = max(0.0, volumeInjCurrCell);

        volumeToInjectList_[cellI] = volumeInjCurrCell;

        while(volumeInjCurrCell > 0.0)
        {
            scalar d = parcelPdf_->sample();
            scalar parcelVolume =
                nPCurrCell
                * mathematicalConstant::pi/6.0*pow3(d);

            // Save the randomly chosen diameter to be assigned later to the
            // parcel
            parcelDiameters_.append(d);
            parcelsInjPerCellList_[cellI]++;
            volumeInjCurrCell -= parcelVolume;
            newVolume += parcelVolume;
        }

        newParcels += parcelsInjPerCellList_[cellI];
    }
}

template<class CloudType>
Foam::scalar Foam::PatchThresholdRndInjection<CloudType>::setNumberOfParticles
(
    const label,
    const scalar,
    const scalar,
    const scalar
)
{
    return nParticlesPerParcel_[currLocalCellI_];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchThresholdRndInjection<CloudType>::PatchThresholdRndInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    patchNames_
    (
        readList<word>(IStringStream(this->coeffDict().lookup("patchNames"))())
    ),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelPdf_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    cellOwners_(),
    referenceField_
    (
        owner.db().objectRegistry::template
        lookupObject<volScalarField>
        (
            this->coeffDict().lookup("referenceField")
        )
    ),
    referenceFieldList_(),
    threshold_(readScalar(this->coeffDict().lookup("thresholdVoidFrac"))),
    volumeToInjectList_(),
    parcelsInjPerCellList_(),
    parcelDiameters_(),
    patchCellsBoxMinCoordinates_(),
    patchCellsBoxDeltaCoordinates_(),
    nParticlesPerParcel_(),
    searchEngine_(this->owner().mesh()),
    currLocalCellI_(0),
    continuousInjection_(this->coeffDict().lookup("allowContinuousInjection")),
    medianParticleVolume_(0.0)
{
    // Find out PDFs median diameter and corresponding volume by sampling the
    // PDF
    scalar dMedian = 0.0;
    label sampleSize = 1000;
    for (int i=0; i < sampleSize; i++)
    {
        dMedian += parcelPdf_->sample();
    }
    dMedian /= sampleSize;
    medianParticleVolume_ = mathematicalConstant::pi/6.0*pow3(dMedian);

    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchThresholdRndInjection<CloudType>::~PatchThresholdRndInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PatchThresholdRndInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::updateMesh()
{
    // Has to be cleared in case the proc has changed
    cellOwners_.clear();

    // Gather cell owners from all patches
    forAll(patchNames_, pI)
    {
        label patchId =
            this->owner().mesh().boundaryMesh().findPatchID(patchNames_[pI]);

        if (patchId < 0)
        {
            FatalErrorIn
            (
                "PatchThresholdRndInjection<CloudType>::"
                "PatchThresholdRndInjection"
                "("
                    "const dictionary&, "
                    "CloudType&"
                ")"
            )   << "Requested patch " << patchNames_[pI] << " not found" << nl
                << "Available patches are: "
                << this->owner().mesh().boundaryMesh().names()
                << nl << exit(FatalError);
        }

        const polyPatch& patch = this->owner().mesh().boundaryMesh()[patchId];

        cellOwners_.append(patch.faceCells());
    }

    cellOwners_.sort();
    label allPatchesSize = cellOwners_.size();
    volumeToInjectList_.resize(allPatchesSize);
    parcelsInjPerCellList_.resize(allPatchesSize);
    patchCellsBoxMinCoordinates_.resize(allPatchesSize);
    patchCellsBoxDeltaCoordinates_.resize(allPatchesSize);
    referenceFieldList_.resize(allPatchesSize);
    nParticlesPerParcel_.resize(allPatchesSize);

    scalar nMaxParticlesPerParcel =
        this->owner().constProps().nMaxParticlesPerParcel();
    scalar nMinParcelsPerCell = this->owner().constProps().nMinParcelsPerCell();

    forAll(cellOwners_, cellI)
    {
        label cellOwner = cellOwners_[cellI];

        labelList cellPointsInd = this->owner().mesh().cellPoints()[cellOwner];
        List<point> cellPoints;

        // Find out the min and max coordinates of all the cell vertices.
        // By doing this we virtually put the cell into a hexahedron
        // (rectangular cuboid)
        vector maxCoord(-vector::max);
        vector minCoord(vector::max);

        forAll(cellPointsInd, i)
        {
            point pt = this->owner().mesh().points()[cellPointsInd[i]];
            cellPoints.append(pt);
        }

        minCoord = min(cellPoints);
        maxCoord = max(cellPoints);

        patchCellsBoxMinCoordinates_[cellI] = minCoord;
        patchCellsBoxDeltaCoordinates_[cellI] = maxCoord - minCoord;

        referenceFieldList_[cellI] = referenceField_[cellOwner];

        // Calculate such nParticle for each parcel that at least a number of
        // nMinParcelsPerCell parcels of the median size fits into each cell.
        scalar cellVolume = this->owner().mesh().cellVolumes()[cellOwner];
        scalar voidFracMedian = medianParticleVolume_ / cellVolume;
        label nP = round(threshold_ / (voidFracMedian*nMinParcelsPerCell));
        if (nP < 1.0)
        {
            nP = 1.0;
        }
        else if (nP > nMaxParticlesPerParcel)
        {
            nP = nMaxParticlesPerParcel;
        }

        nParticlesPerParcel_[cellI] = nP;
    }
}


template<class CloudType>
Foam::scalar Foam::PatchThresholdRndInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    if (cellOwners_.size() > 0)
    {
        // At least one parcel will be injected if this method is called
        while (parcelsInjPerCellList_[currLocalCellI_] == 0)
        {
            currLocalCellI_++;
        }

        cellOwner = cellOwners_[currLocalCellI_];

        rndPositionInPatchCell(currLocalCellI_, cellOwner, position);

        // Update the list with the number of parcels to inject
        parcelsInjPerCellList_[currLocalCellI_]--;
    }
    else
    {
        cellOwner = -1;
        // dummy position
        position = pTraits<vector>::max;
    }
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // particle diameter of the parcel was defined in prepareForNextTimeStep()
    parcel.d() = parcelDiameters_[parcelI];

    label nP = nParticlesPerParcel_[currLocalCellI_];
    nP = this->owner().rndGen().integer(round(0.9*nP), round(1.1*nP));
    parcel.nParticle() = (nP < 1.0) ? 1.0 : nP;

    parcel.lastCell() = parcel.cell();
}


template<class CloudType>
bool Foam::PatchThresholdRndInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::PatchThresholdRndInjection<CloudType>::validInjection(const label)
{
    return true;
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::updateAndInjectParcelIfPatch
(
    typename CloudType::parcelType* origParcel,
    typename CloudType::parcelType::trackData& td,
    const scalar dt,
    const label cellI
)
{
    if (!continuousInjection_ || !isInjectionTime())
    {
        return;
    }

    const label lc = origParcel->lastCell();

    if (lc != cellI)
    {
        const label lastCellLocalI = isInConstVoidFracRegion(lc);
        const scalar curParcelVolume =
            origParcel->nParticle()*origParcel->volume();

        // Parcel comes from an injector patch cell
        if (lastCellLocalI != -1)
        {
            const fvMesh& mesh = td.cloud().mesh();

            // Update void fraction of the last cell
            updateVoidFrac(lc, lastCellLocalI, -curParcelVolume);

            // If void fraction in the last cell is under the threshold
            // inject only one parcel. On average, the threshold is restored.
            if (isUnderThreshold(lastCellLocalI))
            {
                // Create a new parcel in the last cell with random position
                vector pos;
                rndPositionInPatchCell(lastCellLocalI, lc, pos);

                typename CloudType::parcelType* pPtr =
                    new typename CloudType::parcelType
                    (
                        td.cloud(),
                        pos,
                        lc
                    );
                pPtr->d() = parcelPdf_->sample();
                pPtr->nParticle() = nParticlesPerParcel_[lastCellLocalI];

                td.cloud().checkParcelProperties(*pPtr, 0.0, false);
                td.cloud().addParticle(pPtr);

                // Update void fraction due to new parcel
                updateVoidFrac
                (
                    lc,
                    lastCellLocalI,
                    pPtr->nParticle()*pPtr->volume()
                );

                // The new parcel is injected when the original parcel has hit
                // a not shared face of a neigbour of the injector patch cell.
                // Alter new step fraction to be of the original parcel minus
                // the half of step fraction traveled inside the neigbour cell.
                // On average the new parcel is injected (since random) in the
                // middle of the injector patch cell. As a result, most probably
                // the new and the original parcel will be "one cell apart" from
                // each other making injection independent of the lagrangian dt.
                scalar deltaStepFrac = dt/mesh.time().deltaTValue();
                pPtr->stepFraction() =
                    origParcel->stepFraction() - 0.5*deltaStepFrac;

                // Evaluate the new parcel in order for it to have a velocity
                // for the current lagrangian time step
                if (dt > ROOTVSMALL)
                {
                    // switchProcessor from trackData object will be used
                    // in this method when setCellVaues is called on the new
                    // parcel. Make sure that it is set to 0 and will be reset
                    // to its original value after velocity of the new parcel
                    // has been evaluated.
                    bool spOrig = td.switchProcessor;
                    td.switchProcessor = 0;

                    pPtr->setCellValues(td, dt, lc);
                    if (td.cloud().cellValueSourceCorrection())
                    {
                        pPtr->cellValueSourceCorrection(td, dt, lc);
                    }
                    pPtr->calc(td, dt, lc);

                    td.switchProcessor = spOrig;
                }

                // Update counters from trackData
                td.parcelsAddedConti()++;
                td.massAddedConti() += pPtr->nParticle() * pPtr->mass();
            }
        }

        // If the new cell of the current parcel also belongs to the
        // injector patch update its void fraction indifferent whether
        // last cell was in const void fraction patch
        updateVoidFracIfPatch
        (
            cellI,
            curParcelVolume
        );

        origParcel->lastCell() = cellI;
    }

    // If parcel is transfered to different proc, cell might not have changed
    // during the last trackToFace. Update void fraction conditionally on the
    // proc using the current cell because in the next step the parcel will be
    // transfered and the information will be lost.
    if (td.switchProcessor)
    {
        updateVoidFracIfPatch
        (
            cellI,
            -(origParcel->nParticle()*origParcel->volume())
        );
    }

}


template<class CloudType>
label Foam::PatchThresholdRndInjection<CloudType>::isInConstVoidFracRegion
(
    const label cellI
) const
{
    return findSortedIndex<SortableList<label> >(cellOwners_, cellI);
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::rndPositionInPatchCell
(
    const label localCellI,
    const label cellOwner,
    vector& position
)
{
    do
    {
        // Inject position is the sum of the lower conner of the cuboid and
        // a random position inside the cuboid.
        position =
            patchCellsBoxMinCoordinates_[localCellI]
          + cmptMultiply
            (
                this->owner().rndGen().vector01(),
                patchCellsBoxDeltaCoordinates_[localCellI]
            );
    }
    // Since the cuboid is larger than the cell, the point might not lie
    // inside the cell. In this case randomize the position again.
    while (!searchEngine_.pointInCell(position, cellOwner));
    const fvMesh& mesh = this->owner().mesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), position);
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::rndPositionInPatchCell
(
    const label cellOwner,
    vector& position
)
{
    label localCellI =
        findSortedIndex<SortableList<label> >(cellOwners_, cellOwner);

    if (localCellI == -1)
    {
        FatalErrorIn
        (
            "void PatchThresholdRndInjection<CloudType>::rndPositionInPatchCell"
        )
            << "Provided cellOwner is not in the patch list"
            << abort(FatalError);
    }

    rndPositionInPatchCell(localCellI, cellOwner, position);
}


template<class CloudType>
const Foam::pdf& Foam::PatchThresholdRndInjection<CloudType>::parcelPdf() const
{
    // Object is a pointer, dereference it by operator()
    return parcelPdf_();
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::updateVoidFrac
(
    const label cellOwner,
    const label localCellI,
    const scalar volume
)
{
    scalar cellVolume = this->owner().mesh().cellVolumes()[cellOwner];
    // Info<< "Before voidFrac[" << cellOwner << "] = " << referenceFieldList_[localCellI]
    //     << nl
    //     << "volume = " << volume << nl
    //     << "cellVolume = " << cellVolume << nl
    //     << "addVoidFrac = " << volume/cellVolume << nl
    //     << "ListVoidFrac = " << referenceFieldList_ << nl
    //     << endl;
    referenceFieldList_[localCellI] += volume/cellVolume;
    // Info<< "After voidFrac[" << cellOwner << "] = " << referenceFieldList_[localCellI]
    //     << endl;
}


template<class CloudType>
void Foam::PatchThresholdRndInjection<CloudType>::updateVoidFracIfPatch
(
    const label cellOwner,
    const scalar volume
)
{
    label localCellI =
        findSortedIndex<SortableList<label> >(cellOwners_, cellOwner);

    if (localCellI != -1)
    {
        updateVoidFrac(cellOwner, localCellI, volume);
    }
}


template<class CloudType>
bool Foam::PatchThresholdRndInjection<CloudType>::isUnderThreshold
(
    const label localCellI
)
{
    if (referenceFieldList_[localCellI] < threshold_)
    {
        return true;
    }

    return false;
}


template<class CloudType>
bool Foam::PatchThresholdRndInjection<CloudType>::isInjectionTime()
{
    const scalar time = this->owner().db().time().value();
    if (time < this->SOI_ || time > timeEnd())
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
