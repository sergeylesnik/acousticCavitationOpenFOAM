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

\*---------------------------------------------------------------------------*/

#include "StochasticDispersionIncomprRAS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionIncomprRAS<CloudType>::
StochasticDispersionIncomprRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionIncomprRASModel<CloudType>(dict, owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionIncomprRAS<CloudType>::
~StochasticDispersionIncomprRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StochasticDispersionIncomprRAS<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::vector Foam::StochasticDispersionIncomprRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    const scalar cps = 0.16432;

    const volScalarField& k = *this->kPtr_;
    const volScalarField& epsilon = *this->epsilonPtr_;

    const scalar UrelMag = mag(U - Uc - UTurb);

    const scalar tTurbLoc = min
    (
        k[celli]/epsilon[celli],
        cps*pow(k[celli], 1.5)/epsilon[celli]/(UrelMag + SMALL)
    );
    
    // Parcel is perturbed by the turbulence
    // Relax the condition regarding time scales (in comparison to the original
    // implementation). Otherwise parcels close to faces (small dt) fulfill
    // the condition whereas parcels near cell centroid do not.
    if (dt < 100*tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0.0;

            scalar sigma = sqrt(2.0*k[celli]/3.0);
            vector dir = 2.0*this->owner().rndGen().vector01() - vector::one;
            dir /= mag(dir) + SMALL;

            // Numerical Recipes... Ch. 7. Random Numbers...
            scalar x1 = 0.0;
            scalar x2 = 0.0;
            scalar rsq = 10.0;
            while ((rsq > 1.0) || (rsq == 0.0))
            {
                x1 = 2.0*this->owner().rndGen().scalar01() - 1.0;
                x2 = 2.0*this->owner().rndGen().scalar01() - 1.0;
                rsq = x1*x1 + x2*x2;
            }

            scalar fac = sqrt(-2.0*log(rsq)/rsq);

            fac *= mag(x1);

            UTurb = sigma*fac*dir;

        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = vector::zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
