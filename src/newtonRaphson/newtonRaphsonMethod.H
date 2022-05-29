#ifndef newtonRaphsonMethod_H
#define newtonRaphsonMethod_H

#include "fvCFD.H"
#include "blockMUMPSSolver.H"
#include "waveNumTableInterpolator.H"
#include "newtonRaphsonControl.H"

namespace Foam
{

class newtonRaphsonMethod
{

    // Private data

    newtonRaphsonControl& ctrl_;
    autoPtr<blockMUMPSSolver>& linearSolver_;
    const fvMesh& mesh_;
    fvBlockMatrix<vector2>& jacobianEqn_;
    volScalarField& PAcRe_;
    volScalarField& PAcIm_;
    volScalarField& kSqrRe_;
    volScalarField& kSqrIm_;
    waveNumCalc& kSqrInterp_;
    // Temporary field which holds the solution (complex pressure) and suits
    // the instantiation of the fvBlockMatrix. Used as a dummy field,
    // because it is overriden every time the equations are inserted.
    volVector2Field blockPAc_;

public:

    // Constructors

        newtonRaphsonMethod
        (
            newtonRaphsonControl& ctrl,
            autoPtr<blockMUMPSSolver>& directSolver,
            fvBlockMatrix<vector2>& jacobian,
            volScalarField& PAcRe,
            volScalarField& PAcIm,
            volScalarField& kSqrRe,
            volScalarField& kSqrIm,
            waveNumCalc& kSqrInterp
        );

    // Member functions

        void solve
        (
            volScalarField& dPAcRe,
            volScalarField& dPAcIm
        );

        // Build real and imaginary parts of Helmholtz equation and insert them
        // into fvBlockMatrix. The solution field which is referenced at the
        // matrix instantiation will be overriden.
        void buildEqnsAndInsert(fvBlockMatrix<vector2>& blockM);

        scalar scaledL2NormError(fvBlockMatrix<vector2>& matrix);

        scalar weightedL2NormError
        (
            const Field<vector2>& e,
            const Field<vector2>& x
        );
};

} // End namespace Foam


#endif // newtonRaphsonMethod_H
