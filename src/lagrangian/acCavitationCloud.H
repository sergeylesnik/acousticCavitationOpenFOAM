#ifndef ACCAVITATIONCLOUD_H
#define ACCAVITATIONCLOUD_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class acCavitationCloud Declaration
\*---------------------------------------------------------------------------*/

class acCavitationCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        acCavitationCloud(const acCavitationCloud&);

        //- Disallow default bitwise assignment
        void operator=(const acCavitationCloud&);


public:

    //- Runtime type information
    TypeName("acCavitationCloud");

    // Constructors

        //- Null constructor
        acCavitationCloud();


    //- Destructor
    virtual ~acCavitationCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
