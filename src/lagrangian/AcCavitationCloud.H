#ifndef ACCAVITATIONBUBBLECLOUD_H
#define ACCAVITATIONBUBBLECLOUD_H

// Implementation based on ThermoCloudTemplate.H

#include "KinematicCloudTemplate.H"
#include "acCavitationCloud.H"
#include "acCavitationParticleForces.H"
#include "interpolation2DTable.H"
#include "meshSearch.H"

#include "stepSmoother.H"
#include "linearDampingWaveNumber.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pdf;

/*---------------------------------------------------------------------------*\
                        Class AcCavitationCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class AcCavitationCloud
:
    public KinematicCloud<ParcelType>,
    public acCavitationCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        AcCavitationCloud(const AcCavitationCloud&);

        //- Disallow default bitwise assignment
        void operator=(const AcCavitationCloud&);


protected:

    // Protected data

        //- Acoustic cavitation parcel constant properties
        typename ParcelType::constantProperties constProps_;

        // Reference to the carrier liquid fields

            //- Acoustic pressure amplitude
            const volScalarField& PAc_;

            //- Return carrier phase of  acoustic pressure amplitude
            const volScalarField& argPAc_;

            //- Return carrier gradient of acoustic pressure amplitude
            const volVectorField& G_;

            //- Return carrier phase of gradient of acoustic pressure
            //  amplitude
            const volVectorField& argG_;

        // Modelling options

            //- Include damping due to viscous friction
            Switch dampingViscous_;

            //- Include damping due to thermal dissipation
            Switch dampingThermal_;

            //- Include damping due to acoustic radiation
            Switch dampingRadiation_;

            //- Dictionary for initial domain injection
            const dictionary& injectDomainInitDict_;

            //- Switch for inject into the whole domain at the initilization
            Switch injectDomainInit_;

            //- Dict which holds a Switch to identify whether init injection
            //  already happened. Helpful when restarting a run not from the
            //  zero'th iteration.
            IOdictionary injectDomainFlagDict_;




        // Damping. Use volScalarField in order to be viewable in paraview.

            //-
            volScalarField PiViscous_;

            //-
            volScalarField PiThermal_;

            //-
            volScalarField PiRadiation_;

            //-
            volScalarField PiTotal_;

            //-
            volScalarField parcelsPerCell_;

            // Tables

                autoPtr<interpolation2DTable<scalar> > PiViTablePtr_;

                autoPtr<interpolation2DTable<scalar> > PiThTablePtr_;

        //- Basic and acoustic cavitation forces
        acCavitationParticleForces acCavForces_;

        //- Number of parcels injected continuously
        label parcelsAddedConti_;



    // Protected member functions

        void injectParcelInCell
        (
            const fvMesh& mesh,
            const vector& deltaCoord, 
            const vector& minCoord, 
            const meshSearch& searchEngine,
            const autoPtr<pdf>& parcelPdf,
            const label cellI,
            const scalar VCellI,
            scalar nP, 
            label& injCounter,
            scalar& voidFracCellI
        );

        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolveCloud();

            //- Post-evolve
            void postEvolve();


public:

    // Constructors

        //- Construct given carrier gas fields
        AcCavitationCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const volScalarField& PAc,
            const volScalarField& argPAc,
            const volVectorField& G,
            const volVectorField& argG,
            const dimensionedVector& g,
            bool readFields = true
        );


    //- Destructor
    virtual ~AcCavitationCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            //- Return the constant properties
            inline const typename ParcelType::constantProperties&
                constProps() const;

            // Sub-models

            // References to the carrier liquid fields

                //- Return carrier acoustic pressure amplitude
                inline const volScalarField& PAc() const;

                //- Return carrier phase of  acoustic pressure amplitude
                inline const volScalarField& argPAc() const;

                //- Return carrier gradient of acoustic pressure amplitude
                inline const volVectorField& G() const;

                //- Return carrier phase of gradient of acoustic pressure
                //  amplitude
                inline const volVectorField& argG() const;


            // Modelling options

                //- Acoustic radiation damping flag
                inline bool dampingViscous() const;

                //- Acoustic radiation damping flag
                inline bool dampingThermal() const;

                //- Acoustic radiation damping flag
                inline bool dampingRadiation() const;


            // Damping

                //-
                inline DimensionedField<scalar, volMesh>& PiViscous();

                //-
                inline DimensionedField<scalar, volMesh>& PiThermal();

                //-
                inline DimensionedField<scalar, volMesh>& PiRadiation();

                inline tmp<volScalarField> PiTotal(volScalarField& PAc);


            // Forces

                //- Return cavitation particle forces
                inline acCavitationParticleForces& acCavForces();

            //- Parcels injected continuously
            inline label parcelsAddedConti() const;


        // Check

            //- Print cloud information
            void info(typename ParcelType::trackData& td) const;


        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, move)
            void evolve();

            //- Inject into the whole domain up to a specified threshold once at
            //  the initilaization
            void injectDomainInit();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);
            virtual void updateMesh();


        // Acoustic cavitation coupling

            //- Compute imaginary part of the wave number
            void computeImagWaveNumber
            (
                DimensionedField<scalar, volMesh>& kSqrIm,
                DimensionedField<scalar, volMesh>& kSqrRe,
                const DimensionedField<scalar, volMesh>& PAc,
                const DimensionedField<scalar, volMesh>& voidFrac,
                linearDampingWaveNumber& ldwn,
                stepSmoother& ss
            );

            //- Reset all fields regarding acoustics
            void resetAcousticFields();

            //- Use this instead of base nParticlesPerCell to get rid of
            //  polluting Pout
            labelList nParcelsPerCell() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AcCavitationCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AcCavitationCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
