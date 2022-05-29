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
    Foam::PatchThresholdRndInjection

Description
    Patch injection

    - User specifies
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Initial parcel velocity
      - Injection volume flow rate
    - Parcel diameters obtained by PDF model
    - Parcels injected at cell centres adjacent to patch

SourceFiles
    PatchThresholdRndInjection.C

\*---------------------------------------------------------------------------*/

#ifndef PatchThresholdRndInjection_H
#define PatchThresholdRndInjection_H

#include "InjectionModel.H"
#include "meshSearch.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class DataEntry;

class pdf;

/*---------------------------------------------------------------------------*\
                       Class PatchThresholdRndInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchThresholdRndInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Names of patches
        const List<word> patchNames_;

        //- Injection duration [s]
        const scalar duration_;

        //- Parcel size PDF model
        const autoPtr<pdf> parcelPdf_;

        //- Cell owners
        SortableList<label> cellOwners_;

        //- Reference Field
        const volScalarField& referenceField_;

        //- List with reference Field values belonging to the patch
        scalarList referenceFieldList_;

        //- Threshold
        const scalar threshold_;

        //- Volume to inject per time step
        scalarList volumeToInjectList_;

        labelList parcelsInjPerCellList_;
        DynamicList<scalar> parcelDiameters_;
        List<vector> patchCellsBoxMinCoordinates_;
        List<vector> patchCellsBoxDeltaCoordinates_;
        labelList nParticlesPerParcel_;
        meshSearch searchEngine_;

        //- Current cell index from the local patch cell list
        label currLocalCellI_;

        Switch continuousInjection_;

        scalar medianParticleVolume_;


protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step relative to SOI
        label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;

        //- Volume of parcels to introduce over the time step relative to SOI
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


        void prepareForNextTimeStep
        (
            const scalar time,
            label& newParcels,
            scalar& newVolume
        );

        //- Set number of particles to inject given parcel properties
        virtual scalar setNumberOfParticles
        (
            const label,
            const scalar,
            const scalar,
            const scalar
        );

public:

    //- Runtime type information
    TypeName("PatchThresholdRndInjection");


    // Constructors

        //- Construct from dictionary
        PatchThresholdRndInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~PatchThresholdRndInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);

        //- Inject a single parcel during the lagrangian dt and update void frac
        void updateAndInjectParcelIfPatch
        (
            typename CloudType::parcelType* origParcel,
            typename CloudType::parcelType::trackData& td,
            const scalar dt,
            const label cellI
        );

        //- Check if a cell is in the const void fraction region and return the
        //  corresponding index from the local patch cell list
        label isInConstVoidFracRegion(const label cellI) const;

        void rndPositionInPatchCell
        (
            const label localCellI,
            const label cellOwner,
            vector& position
        );

        void rndPositionInPatchCell(const label cellOwner, vector& position);

        inline const pdf& parcelPdf() const;

        void updateVoidFrac
        (
            const label cellOwner,
            const label localCellI,
            const scalar volume
        );

        void updateVoidFracIfPatch(const label cellOwner, const scalar volume);

        bool isUnderThreshold(const label localCellI);

        bool isInjectionTime();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatchThresholdRndInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
