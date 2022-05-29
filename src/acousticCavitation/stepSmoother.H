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

#ifndef stepSmoother_H
#define stepSmoother_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class stepSmoother
{

protected:

    // Protected Data

    const scalar spread_;

    scalar lowerLimit_;

    scalar upperLimit_;

public:

    // Constructor

        stepSmoother(const scalar spread)
        :
            spread_(spread),
            lowerLimit_(0.0),
            upperLimit_(0.0)
        {}

    // Member Functions

        void setLimits(scalar midValue)
        {
            lowerLimit_ = (1 - spread_)*midValue;
            upperLimit_ = (1 + spread_)*midValue;
        }

        scalar clamp(scalar x, scalar lowerlimit, scalar upperlimit)
        {
            if (x < lowerlimit)
            x = lowerlimit;
            if (x > upperlimit)
            x = upperlimit;
            return x;
        }

        scalar smoothStep(scalar edge0, scalar edge1, scalar x)
        {
            // Scale, bias and saturate x to 0..1 range
            x = clamp((x - edge0)/(edge1 - edge0), 0.0, 1.0);
            // Evaluate polynomial
            return x*x*(3 - 2*x);
        }

        scalar smoothStep(scalar midValue, scalar x)
        {
            setLimits(midValue);
            return smoothStep(lowerLimit_, upperLimit_, x);
        }

        scalar smoothStep(scalar x)
        {
            return smoothStep(lowerLimit_, upperLimit_, x);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
