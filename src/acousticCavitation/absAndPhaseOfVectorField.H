#ifndef absAndPhaseOfVectorField_H
#define absAndPhaseOfVectorField_H

namespace Foam
{

void absAndPhaseOfVectorField
(
    volVectorField& fieldAbs,
    volVectorField& fieldPhase,
    const volVectorField &fieldReal,
    const volVectorField &fieldImag
)
{

    // sqr() and sqrt() work only for volScalarField (no element-wise
    // methods for volVectorField). Thus, loop over components.
    for (int dir = 0; dir < 3; dir++)
    {
        volScalarField cmptReal(fieldReal.component(dir));
        volScalarField cmptImag(fieldImag.component(dir));
        fieldAbs.replace
        (
            dir,
            sqrt(sqr(cmptReal) + sqr(cmptImag))
        );

        // atan2() works properly only for Field<> (for GeometricField<>
        // it delivers Field<>).
        // Thus, replace internal and boundary fields separately.
        fieldPhase.internalField().replace(dir, atan2(cmptImag, cmptReal));
        fieldPhase.boundaryField().replace
        (
            dir,
            atan2(cmptImag.boundaryField(), cmptReal.boundaryField())
        );
    }
    // In case of wedge, constrain the circumferential direction of the fields.
    const fvMesh& mesh = fieldAbs.mesh();
    meshTools::constrainDirection(mesh, mesh.geometricD(), fieldAbs);
    meshTools::constrainDirection(mesh, mesh.geometricD(), fieldPhase);
    fieldAbs.correctBoundaryConditions();
    fieldPhase.correctBoundaryConditions();
}


} // End namespace Foam

#endif // absAndPhaseOfVectorField_H
