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
    Foam::StochasticDispersionIncomprRAS

Description
    The velocity is perturbed in random direction, with a
    Gaussian random number distribution with variance sigma.
    Source:
    [1] A. D. Gosman and E. Ioannides, “Aspects of Computer Simulation of
    Liquid-Fueled Combustors,” J. Energy, vol. 7, no. 6, pp. 482–490, Nov.
    1983, doi: 10.2514/3.62687.

\*---------------------------------------------------------------------------*/

#ifndef StochasticDispersionIncomprRAS_H
#define StochasticDispersionIncomprRAS_H

#include "DispersionIncomprRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class StochasticDispersionIncomprRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StochasticDispersionIncomprRAS
:
    public DispersionIncomprRASModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("StochasticDispersionIncomprRAS");


    // Constructors

        //- Construct from components
        StochasticDispersionIncomprRAS
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~StochasticDispersionIncomprRAS();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        virtual bool active() const;

        //- Update (disperse particles)
        virtual vector update
        (
            const scalar dt,
            const label celli,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StochasticDispersionIncomprRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
