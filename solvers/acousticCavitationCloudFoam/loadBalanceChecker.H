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

#ifndef loadBalanceChecker_H
#define loadBalanceChecker_H

#include "AcCavitationCloud.H"
#include "acCavitationParcel.H"
#include "labelIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class loadBalanceChecker
{
    // Private data

        const fvMesh& mesh_;

        const AcCavitationCloud<acCavitationParcel>& cloud_;

        IOdictionary dynamicMeshDict_;

        const dictionary& loadBalanceCoeffs_;

        scalar imbalanceTrigger_;

        const dictionary& parMetisByParcelsCoeffs_;

        word weightsFileName_;

        labelIOField cellIOWeights_;
        
        scalar injectedParcelsWeightAmplifier_;

public:

    // Constructors

        loadBalanceChecker
        (
            const fvMesh& mesh,
            const AcCavitationCloud<acCavitationParcel>& cloud
        );

    // Member Funcitons

        bool isImba();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //