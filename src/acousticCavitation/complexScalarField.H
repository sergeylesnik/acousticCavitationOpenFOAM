#ifndef complexScalarField_H
#define complexScalarField_H

#include "fvCFD.H"

class complexScalarField
{
    // Private Data
    volScalarField& Re_;
    volScalarField& Im_;
    volScalarField& Abs_;

public:

    // Constructors
    complexScalarField
    (
        volScalarField& Re,
        volScalarField& Im,
        volScalarField& Abs
    )
    :
    Re_(Re),
    Im_(Im),
    Abs_(Abs)
    {}

    // Public Methods
    volScalarField& Re()
    {
        return Re_;
    }

    volScalarField& Im()
    {
        return Im_;
    }

    volScalarField& Re() const
    {
        return Re_;
    }

    volScalarField& Im() const
    {
        return Im_;
    }

    volScalarField& computeAbsolute()
    {
        Abs_ = sqrt(pow(Re_,2) + pow(Im_,2));
        return Abs_;
    }
};

#endif // complexScalarField_H

