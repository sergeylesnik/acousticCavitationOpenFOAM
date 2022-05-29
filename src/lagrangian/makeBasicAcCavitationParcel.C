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

#include "acCavitationParcel.H"

// Kinematic
// #include "makeParcelDispersionModels.H"
#include "makeParcelDragModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelPostProcessingModels.H"

#include "MagnaudetDrag.H"

// Acoustic cavitation
#include "PatchThresholdRndInjection.H"

// Dispersiond model for incompressible turbulence RAS
#include "NoDispersion.H"
#include "StochasticDispersionIncomprRAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    // makeParcelDispersionModels(acCavitationParcel);
    makeParcelDragModels(acCavitationParcel);
    makeParcelInjectionModels(acCavitationParcel);
    makeParcelPatchInteractionModels(acCavitationParcel);
    makeParcelPostProcessingModels(acCavitationParcel);

    // Additional acoustic cavitation sub-models
    makeInjectionModelType                                                    \
    (                                                                         \
        PatchThresholdRndInjection,                                           \
        KinematicCloud,                                                       \
        acCavitationParcel                                                    \
    )
    
    makeDragModelType(MagnaudetDrag, KinematicCloud, acCavitationParcel);

    // Dispersion model. The only change: compressible -> incompressible
    makeDispersionModel(KinematicCloud<acCavitationParcel>);

    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        DispersionIncomprRASModel<KinematicCloud<acCavitationParcel> >,       \
        0                                                                     \
    );

    makeDispersionModelType                                                   \
    (                                                                         \
        NoDispersion,                                                         \
        KinematicCloud,                                                       \
        acCavitationParcel                                                    \
    );

    makeDispersionModelType                                                   \
    (                                                                         \
        StochasticDispersionIncomprRAS,                                       \
        KinematicCloud,                                                       \
        acCavitationParcel                                                    \
    );
};


// ************************************************************************* //
