#include "acCavitationParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(acCavitationParcel, 0);
    defineParticleTypeNameAndDebug(acCavitationParcel, 0);
    defineParcelTypeNameAndDebug(acCavitationParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::acCavitationParcel::acCavitationParcel
(
    const acCavitationParcel& p
)
:
    KinematicParcel<acCavitationParcel>(p),
    lastCell_(p.cell()),
    nSmallDeltaStepFractions_(p.nSmallDeltaStepFractions())
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::acCavitationParcel::~acCavitationParcel()
{}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

const Foam::vector Foam::acCavitationParcel::calcVelocity
(
    trackData& td,
    const scalar dt,
    const label cellI,
    const scalar mu,
    const scalar dAv,
    const vector& U,
    const scalar rho,
    const scalar mass,
    const vector& Su,
    vector& dUTrans
) const
{
    const polyMesh& mesh = this->cloud().pMesh();

    // scalar utc;
    // scalar As = this->areaS(dAv);
    IntegrationScheme<vector>::integrationResult Ures;
    vector Unew = U;

    // Momentum source due to particle forces
    scalar massEff = mass + td.cloud().acCavForces().massAdd(rhoc_, VAv_);

    const vector FCoupled =
        mass
      * td.cloud().acCavForces().calcCoupled(cellI, dt, rhoc_, rho, Uc_, U);
    const vector FNonCoupled =
        mass
      * td.cloud().acCavForces().calcNonCoupled
        (
            cellI,
            dt,
            rhoc_,
            rho,
            Uc_,
            U,
            FBjPri_/mass
        );

    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Sergey Lesnik Bugfix 13.07.2020
    // Avoid zero drag coefficient (utc) for zero Reynolds.
    scalar Re = max(SMALL, this->Re(U, dAv, rhoc_, muc_));
    scalar RePrev;
    scalar SpDrag;

    // Drag explicit source Sp(Fd) might be off several of its multiplies.
    // Loop until bubble Re and thus Cd, Sp(Fd) converge.
    do
    {
        // Momentum transfer coefficient
        // utc = td.cloud().drag().utc(Re, dAv, mu) + ROOTVSMALL;
        scalar Cd = td.cloud().drag().Cd(Re);
        SpDrag = Cd*Re*mu*mathematicalConstant::pi*dAv/8.0 + ROOTVSMALL;

        // Update velocity - treat as 3-D
        // const vector apOld =
        //     Uc_ + (FCoupled + FNonCoupled + Su)/(utc*As)*mass/massEff;
        // const scalar bpOld = 6.0*utc/(rho*dAv);
        const vector ap = Uc_ + (FCoupled + FNonCoupled + Su)/SpDrag;
        const scalar bp = SpDrag/massEff;

        Ures = td.cloud().UIntegrator().integrate(U, dt, ap, bp);
        Unew = Ures.value();

        RePrev = Re;
        Re = max(SMALL, this->Re(Unew, dAv, rhoc_, muc_));
// if (debug)
// {
// Info<< "Unew = " << Unew << nl
//     << "FNonCoup = " << FNonCoupled << nl
//     << "abs = " << mag(RePrev - Re) << nl
//     << "ap*bp = " << ap*bp << nl
//     // << "apOld*bpOld = " << apOld*bpOld << nl
//     << "bp = " << bp << nl
//     // << "bpOld = " << bpOld << nl
//     << "mass = " << mass << nl
//     << "massEff = " << massEff << nl
//     << "mAv_ = " << rho_*VAv_ << nl
//     << "VAv_ = " << VAv_ << nl
//     << "Vp = " << volume(d_) << nl
//     << "V(dAv) = " << volume(dAv) << nl
//     << "SpDrag = " << SpDrag << nl
//     << "utc*As = " << utc*As << nl << endl;
// }
    }
    while (mag(RePrev - Re) > td.ReAbsTol());

    dUTrans += dt*(SpDrag*(Ures.average() - Uc_) - FCoupled);

    // Apply correction to velocity and dUTrans for reduced-D cases
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::acCavitationParcel::setCellValues
(
    trackData &td,
    const scalar dt,
    const label cellI
)
{
    const fvMesh& mesh = td.cloud().mesh();

    // Shift towards the cell center parcels with very small delta step
    // fraction (happens e.g. at the axis of wedge geometry when parcel hits
    // both wedge patches alternately). Exclude cases where its the last step
    // (stepFrac = 1.0) and shift only if this behaviour occurs some n times
    // for the same parcel.
    scalar deltaStepFraction = dt/mesh.time().deltaTValue();
    label nLimitSmallStepFrac = 10;
    if (deltaStepFraction < 1e-3 && stepFraction() != 1.0)
    {
        nSmallDeltaStepFractions_++;

        if (nSmallDeltaStepFractions_ > nLimitSmallStepFrac)
        {
            position() += 1.0e-1*(mesh.cellCentres()[cellI] - position());

            td.nParcelShifts()++;

            if (nSmallDeltaStepFractions_ == nLimitSmallStepFrac + 1)
            {
                td.nShiftedParcels()++;
            }
        }
    }
    // Reset the counter if the parcel is evaluated the last time
    else if (stepFraction() == 1.0)
    {
        nSmallDeltaStepFractions_ = 0;
    }

    KinematicParcel<acCavitationParcel>::setCellValues(td, dt, cellI);

    PAc_ = td.PAcInterp().interpolate(this->position(), cellI);
    dAv_ = 2 * td.RAvTable()(REqu(), PAc_);
    VAv_ = td.VAvTableInterpolate(REqu(), PAc_);

    // Compute parcel acceleration due to primary Bjerknes force
    scalar argPAc = td.argPAcInterp().interpolate(this->position(), cellI);
    vector G = td.GInterp().interpolate(this->position(), cellI);
    vector argG = td.argGInterp().interpolate(this->position(), cellI);
    scalar Ic = td.IcTableInterpolate(REqu(), PAc_);
    scalar Is = td.IsTableInterpolate(REqu(), PAc_);

    // cos() and sin() do not support Foam::vector. Iterate over components.
    forAll(G, cmptI)
    {
        FBjPri_[cmptI] =
            G[cmptI]
          * (
                Ic*cos(argPAc-argG[cmptI])
              + Is*sin(argPAc-argG[cmptI])
            );
    }

    // Constrain in case of 1D or 2D
    meshTools::constrainDirection(mesh, mesh.solutionD(), FBjPri_);

    if (debug)
    {
        label lc = lastCell();
        Pout<< "switchProc = " << td.switchProcessor << nl
            << "lastCell = " << lc << nl
            << "curCell = " << cellI << nl
            << "stepFraction = " << stepFraction()<< nl
            << "dt = " << dt << nl
            << "deltaTValue = " << mesh.time().deltaTValue() << nl
            << "deltaStepFrac = " << dt/mesh.time().deltaTValue() << nl
            << "origId = " << origId() << nl
            << "origProc = " << origProc() << nl
            << "pos = " << position() << nl
            << "U = " << U() << nl
            << "UTurb = " << UTurb() << nl
            << "tTurb = " << tTurb() << nl
            << endl;
    }

    // Inject a new parcel if the parcel comes from the cell marked as const
    // voidFrac
    if (td.patchRndInjectorPtr())
    {
        td.patchRndInjectorPtr()->
            updateAndInjectParcelIfPatch(this, td, dt, cellI);
    }
}


// calc might be called several times per time step, namely every time parcel
// leaves a cell. dt is the corresponding fraction of the lagrangian time step.
void Foam::acCavitationParcel::calc
(
    trackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar dAv = dAv_;
    const vector U0 = U_;
    const scalar rho0 = rho_;
//    const scalar mass0 = mass();
    // Mass averaged over one acoustic time period.
    const scalar massAv = rho_*VAv_;


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    vector U1 =
        calcVelocity(td, dt, cellI, muc_, dAv, U0, rho0, massAv, Su, dUTrans);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    U_ = U1;
}

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "acCavitationParcelIO.C"

// ************************************************************************* //
