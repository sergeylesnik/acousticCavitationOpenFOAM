#ifndef acCavitationParcel_H
#define acCavitationParcel_H

// Implementation based on basicKinemticParcel.H

#include "KinematicParcel.H"
#include "AcCavitationCloud.H"
#include "interpolation2DTable.H"
#include "PatchThresholdRndInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class acCavitationParcel;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const acCavitationParcel&
);

/*---------------------------------------------------------------------------*\
                   Class acCavitationParcel Declaration
\*---------------------------------------------------------------------------*/

class acCavitationParcel
:
    public KinematicParcel<acCavitationParcel>
{
public:

    //- Class to hold acoustic cavitation bubble constant properties
    class constantProperties
    :
        public KinematicParcel<acCavitationParcel>::constantProperties
    {
        // Private Data

            const dimensionedScalar freq_;

            const dimensionedScalar pInf_;

            const dimensionedScalar sigma_;

            const dimensionedScalar DThermal_;

            const dimensionedScalar gamma_;

            const scalar spread_;

            const scalar nMaxParticlesPerParcel_;

            const scalar nMinParcelsPerCell_;

        public:

            //- Constructor
            constantProperties(const dictionary& parentDict);

            //- Access

                inline const dimensionedScalar& sigma() const;

                inline const dimensionedScalar& pInf() const;

                inline scalar nMaxParticlesPerParcel() const;

                inline scalar nMinParcelsPerCell() const;
    };


    //- Class used to pass Bjerknes forces to the trackToFace function
    class trackData
    :
        public KinematicParcel<acCavitationParcel>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            AcCavitationCloud<acCavitationParcel>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Absolute tolerance for Re of a parcel
            const scalar ReAbsTol_;

            // Interpolators for continuous phase fields

                //- Acoustic pressure amplitude
                const interpolation<scalar>& PAcInterp_;

                //- Phase of acoustic pressure amplitude [rad]
                const interpolation<scalar>& argPAcInterp_;

                //- Gradient of the acoustic pressure amplitude
                const interpolation<vector>& GInterp_;

                //- Phase of the gradient of the acoustic pressure amplitude
                const interpolation<vector>& argGInterp_;


            // Interpolators for 2D tables with cavitation bubble properties

                //-
                const interpolation2DTable<scalar> RAvTable_;

                const interpolation2DTable<scalar> VAvTable_;

                const interpolation2DTable<scalar> IcTable_;

                const interpolation2DTable<scalar> IsTable_;

            //- Pointer to patch random injector
            PatchThresholdRndInjection<KinematicCloud<acCavitationParcel> >*
                patchRndInjectorPtr_;

            //- Parcel counter for the continuous injection
            label parcelsAddedConti_;

            //- Parcel counter for the continuous injection
            scalar massAddedConti_;

            // Statistics for shifting of parcels towards the cell center
            // when the delta step fraction is very small

                //- Number of parcels which were shifted at least once
                label nShiftedParcels_;

                //- Number of total shifts for all parcels
                label nParcelShifts_;

    public:

        // Constructors

            //- Construct from components
            inline trackData
            (
                AcCavitationCloud<acCavitationParcel>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const interpolation<scalar>& PAcInterp,
                const interpolation<scalar>& argPAcInterp,
                const interpolation<vector>& GInterp,
                const interpolation<vector>& argGInterp,
                const vector& g
            );

            //- Construct as a copy
            inline trackData(trackData& td);

        // Member functions

            //- Return access to the owner cloud
            inline AcCavitationCloud<acCavitationParcel>& cloud();

            //- Return const access to the owner cloud
            inline const AcCavitationCloud<acCavitationParcel>& cloud() const;

            //- Return const access to the owner cloud
            inline const constantProperties& constProps() const;

            //- Return const access to the Reynolds absolute tolerance
            inline scalar ReAbsTol() const;

            //- Return const access to the interpolator for continuous
            //  phase acoustic pressure amplitude
            inline const interpolation<scalar>& PAcInterp() const;

            inline const interpolation<scalar>& argPAcInterp() const;

            inline const interpolation<vector>& GInterp() const;

            inline const interpolation<vector>& argGInterp() const;

            //- Const access to the interpolation table of averaged radius
            inline const interpolation2DTable<scalar>& RAvTable() const;

            //- Const access to the interpolation table of averaged volume
            inline const interpolation2DTable<scalar>& VAvTable() const;

            //- Const access to the interpolation table of Ic
            inline const interpolation2DTable<scalar>& IcTable() const;

            //- Const access to the interpolation table of Is
            inline const interpolation2DTable<scalar>& IsTable() const;

            //- Return volume averaged over an acoustic time period interpolated
            //  in the characteristic map of PAc and REqu
            inline scalar VAvTableInterpolate(scalar REqu, scalar PAc) const;

            //- Return V*cos(omegat*t) averaged over an acoustic time period
            //  interpolated in the characteristic map of PAc and REqu
            inline scalar IcTableInterpolate(scalar REqu, scalar PAc) const;

            //- Return V*sin(omegat*t) averaged over an acoustic time period
            //  interpolated in the characteristic map of PAc and REqu
            inline scalar IsTableInterpolate(scalar REqu, scalar PAc) const;

            inline PatchThresholdRndInjection
                <
                    KinematicCloud<acCavitationParcel>
                >*
            patchRndInjectorPtr();

            inline label& parcelsAddedConti();

            inline scalar& massAddedConti();

            inline label& nShiftedParcels();

            inline label& nParcelShifts();
    };


protected:

    // Protected data

        // Parcel properties

            //- Diameter averaged over one time period [m]
            scalar dAv_;

            //- Volume averaged over one acoustic time period [m^3]
            scalar VAv_;

            //- Last cell the parcel was in, needed for patch injection
            label lastCell_;

            //- Number of times the parcel experiences small delta step fraction
            label nSmallDeltaStepFractions_;

        // Cashed cell-based quantities are updated from the continious phase
        // for every parcel.

            //- Primary Bjerknes force [N]
            vector FBjPri_;

            //- Acoustic pressure amplitude [Pa]
            scalar PAc_;


    // Protected member functions

        //- Calculate new particle velocity
        const vector calcVelocity
        (
            trackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar mu,           // local carrier viscosity
            const scalar d,            // diameter
            const vector& U,           // velocity
            const scalar rho,          // density
            const scalar mass,         // mass
            const vector& Su,          // explicit particle momentum source
            vector& dUTrans            // momentum transfer to carrier
        ) const;


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Run-time type information
        TypeName("acCavitationParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Used in the PatchThresholdRndInjection and as for domain initial
        //  injection.
        inline acCavitationParcel
        (
            KinematicCloud<acCavitationParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        inline acCavitationParcel
        (
            KinematicCloud<acCavitationParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        acCavitationParcel
        (
            const Cloud<acCavitationParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        acCavitationParcel(const acCavitationParcel& p);

        //- Construct and return a clone
        autoPtr<acCavitationParcel> clone() const
        {
            return
                autoPtr<acCavitationParcel>
                (
                    new acCavitationParcel(*this)
                );
        }


    //- Destructor
    virtual ~acCavitationParcel();

    // Member functions

        //- Interpolate needed quantities from eulerian to lagrangian
        void setCellValues
        (
            trackData& td,
            const scalar dt,
            const label cellI
        );

        //- Main function for moving the bubble
        void calc
        (
            trackData& td,
            const scalar dt,
            const label cellI
        );

        // Access

            //- Return cavitation bubble equilibrium radius
            scalar REqu() const;

            //- Return cavitation bubble diameter averaged over one time period
            scalar dAv() const;

            //- Return cell-based value of the primary Bjerknes force
            vector FBjPri() const;

            //- Return cell-based value of the acoustic pressure amplitude
            scalar PAc() const;

            //- Return index of the last cell the bubble was in
            inline label lastCell() const;

            //- Return number of small delta step fractions experienced by the
            //  bubble
            inline label nSmallDeltaStepFractions() const;


        // Edit

            //- Return access to index of the last cell the bubble was in
            inline label& lastCell();


        // I-O

            //- Read
            static void readFields(Cloud<acCavitationParcel>& c);

            //- Write
            static void writeFields(const Cloud<acCavitationParcel>& c);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const acCavitationParcel&
        );
};


template<>
inline bool contiguous<acCavitationParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "acCavitationParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
