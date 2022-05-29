#ifndef residualNorm_H
#define residualNorm_H

#include "runTimeSelectionTables.H"
#include "fvc.H"
#include "fvBlockMatrix.H"

namespace Foam
{

class residualNorm
{
    public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            residualNorm,
            dict,
            (),
            ()
        );

    TypeName("residualNorm");

        static autoPtr<residualNorm> New(const dictionary& dict);

        virtual scalar residual(const fvBlockMatrix<vector2>& matrix) = 0;

        void cmpxResidualField
        (
            const fvBlockMatrix<vector2>& matrix,
            vector2Field& res
        );

        virtual ~residualNorm(){}

};

class scaledL1Norm
:
public residualNorm
{
    public:

        TypeName("scaledL1Norm");

        scalar residual(const fvBlockMatrix<vector2>& matrix);
};

class normalizedL1Norm
:
public residualNorm
{
    public:

        TypeName("normalizedL1Norm");

        scalar residual(const fvBlockMatrix<vector2>& matrix);
};

class scaledL2Norm
:
public residualNorm
{
    public:

        TypeName("scaledL2Norm");

        scalar residual(const fvBlockMatrix<vector2>& matrix);
};

class normalizedL2Norm
:
public residualNorm
{
    public:

        TypeName("normalizedL2Norm");

        scalar residual(const fvBlockMatrix<vector2>& matrix);
};



} // End namespace Foam

#endif // residualNorm_H
