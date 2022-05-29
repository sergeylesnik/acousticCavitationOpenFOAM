#include "bugFixedPimpleControl.H"
#include "fvc.H"

void Foam::bugFixedPimpleControl::calcTransientConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& rAU,
    const fvVectorMatrix& ddtUEqn
) const
{
    // Store necessary data for this velocity field
    const word& UName = U.name();

    // Check whether the fields are present in the list
    if (!indices_.found(UName))
    {
        // Get current index as size of the indices list (before insertion)
        const label i = indices_.size();

        // Insert the index into the hash table
        indices_.insert(UName, i);

        // Extend lists
        aCoeffPtrs_.resize(indices_.size());
        faceUPtrs_.resize(indices_.size());

        // Double check whether the fields have been already set
        if (!aCoeffPtrs_.set(i) && !faceUPtrs_.set(i))
        {
            aCoeffPtrs_.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "aCoeff." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );

            faceUPtrs_.set
            (
                i,
                new surfaceVectorField
                (
                    IOobject
                    (
                        "faceU." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("zero", dimVelocity, vector::zero)
                )
            );

            // Bug fix, Sergey Lesnik, 19.02.2021: Initialize faceU.oldTime
            // using U.oldTime. Needed if oldTime flow field differs from
            // the current one (e.g. restart).
            faceUPtrs_[i].oldTime() = fvc::interpolate(U.oldTime());
        }
        else if (!aCoeffPtrs_.set(i) || !faceUPtrs_.set(i))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Either aCoeff or faceU is allocated for field " << UName
                << " while the other is not." << nl
                << " This must not happen in transient simulation. Make sure"
                << " that functions aiding consistency are called in the right"
                << " order (first flux and then velocity reconstruction)."
                << exit(FatalError);
        }
    }
    else
    {
        // Index has been set for this field, so the fields must be there as
        // well. Check and report an error if they are not allocated
        if (!aCoeffPtrs_.set(indices_[UName]))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Index is set, but the aCoeff field is not allocated for "
                << UName << "." << nl
                << "This should not happen for transient simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
        else if (!faceUPtrs_.set(indices_[UName]))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Index is set, but the faceU field is not allocated for "
                << UName << "." << nl
                << "This should not happen for transient simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
    }

    // Algorithm:
    // 1. Update flux and aCoeff due to ddt discretisation
    // 2. Update flux and aCoeff due to under-relaxation
    // 3. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get index from the hash table
    const label i = indices_[UName];

    // Get fields that will be updated
    surfaceScalarField& aCoeff = aCoeffPtrs_[i];
    surfaceVectorField& faceU = faceUPtrs_[i];

    // Update face interpolated velocity field. Note: handling of oldTime faceU
    // fields happens in ddt scheme when calling ddtConsistentPhiCorr inside
    // addDdtFluxContribution
    faceU = fvc::interpolate(U);

// Info<< "oldTime = " << faceU.oldTime() << nl
//     << "faceU = " << faceU << endl;
// faceU.oldTime() = fvc::interpolate(U.oldTime());
// Info<< "oldTime = " << faceU.oldTime() << endl;

    // Interpolate original rAU on the faces
    const surfaceScalarField rAUf = fvc::interpolate(rAU);

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (faceU & mesh_.Sf());

    // Initialize aCoeff to 1
    aCoeff = dimensionedScalar("one", dimless, 1.0);

    // STAGE 1: consistent ddt discretisation handling
    addDdtFluxContribution(phi, aCoeff, faceU, U, rAUf, ddtUEqn);

    // STAGE 2: consistent under-relaxation handling
    addUnderRelaxationFluxContribution(phi, aCoeff, U);

    // STAGE 3: scale the flux and correct it at the boundaries
    phi /= aCoeff;
    correctBoundaryFlux(phi, U);
}
