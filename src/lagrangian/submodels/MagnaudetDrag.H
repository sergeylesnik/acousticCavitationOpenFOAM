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
    Foam::MagnaudetDrag

Description
    Drag model based on assumption of oscillating bubble with high frequency.
    Source:
    Magnaudet, J. & Legendre, D. The viscous drag force on a spherical bubble
    with a time-dependent radius. Phys. Fluids 10, 550–554 (1998).

\*---------------------------------------------------------------------------*/

#ifndef MagnaudetDrag_H
#define MagnaudetDrag_H

#include "DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class MagnaudetDrag Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MagnaudetDrag
:
    public DragModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("MagnaudetDrag");


    // Constructors

        //- Construct from dictionary
        MagnaudetDrag
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~MagnaudetDrag();


    // Member Functions

        //- Flag to indicate whether model activates drag model
        bool active() const;

        //- Return drag coefficient
        scalar Cd(const scalar Re) const;

        //- Return drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;

        //- Return momentum transfer coefficient
        //  Drag force per unit particle surface area = utc(U - Up)
        scalar utc(const scalar Re, const scalar d, const scalar mu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MagnaudetDrag.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
