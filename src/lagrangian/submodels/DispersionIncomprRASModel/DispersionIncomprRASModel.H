/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DispersionIncomprRASModel

Description

\*---------------------------------------------------------------------------*/

#ifndef DispersionIncomprRASModel_H
#define DispersionIncomprRASModel_H

#include "DispersionModel.H"
#include "RASModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DispersionIncomprRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DispersionIncomprRASModel
:
    public DispersionModel<CloudType>
{
protected:

    // Protected data

        //- Reference to the incompressible turbulence model
        const incompressible::RASModel& turbulence_;

        // Locally cached turbulence fields

            //- Turbulence k
            const volScalarField* kPtr_;

            //- Take ownership of the k field
            bool ownK_;

            //- Turbulence epsilon
            const volScalarField* epsilonPtr_;

            //- Take ownership of the epsilon field
            bool ownEpsilon_;


public:

    //- Runtime type information
    TypeName("DispersionIncomprRASModel");


    // Constructors

        //- Construct from components
        DispersionIncomprRASModel
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~DispersionIncomprRASModel();


    // Member Functions

        //- Cache carrier fields
        virtual void cacheFields(const bool store);

        //- Return const access to the turbulence model
        const incompressible::RASModel& turbulence() const
        {
            return turbulence_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DispersionIncomprRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
