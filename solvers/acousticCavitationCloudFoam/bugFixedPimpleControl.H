#ifndef bugFixedPimpleControl_H
#define bugFixedPimpleControl_H

#include "fvMesh.H"
#include "pimpleControl.H"

namespace Foam
{

class bugFixedPimpleControl
:
    public pimpleControl
{

public:

    bugFixedPimpleControl(fvMesh& mesh)
    :
    pimpleControl(mesh)
    {}

    ~bugFixedPimpleControl()
    {}

    // Override the base method with bug fix
    void calcTransientConsistentFlux
    (
        surfaceScalarField& phi,
        const volVectorField& U,
        const volScalarField& rAU,
        const fvVectorMatrix& ddtUEqn
    ) const;

};

} // End namespace Foam

#endif
