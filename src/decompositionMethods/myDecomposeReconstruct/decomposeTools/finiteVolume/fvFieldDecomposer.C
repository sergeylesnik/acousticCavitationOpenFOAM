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

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const label sizeBeforeMapping,
    const unallocLabelList& addressingSlice,
    const label addressingOffset
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    directAddressing_(addressingSlice)
{
    forAll (directAddressing_, i)
    {
        // Subtract one to align addressing.  HJ, 5/Dec/2001
        directAddressing_[i] -= addressingOffset + 1;
    }
}


fvFieldDecomposer::processorVolPatchFieldDecomposer::
processorVolPatchFieldDecomposer
(
    const fvMesh& mesh,
    const unallocLabelList& addressingSlice
)
:
    sizeBeforeMapping_(mesh.nCells()),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    // Cannot use weights due to parallel comms in load balancing
    // Use 0.5 instead, as this is not precise mapping
    // HJ, 21/Oct/2018

    const labelList& own = mesh.faceOwner();
    const labelList& neighb = mesh.faceNeighbour();

    forAll (addressing_, i)
    {
        // Subtract one to align addressing.  HJ, 5/Dec/2001
        label ai = mag(addressingSlice[i]) - 1;

        if (ai < neighb.size())
        {
            // This is a regular face. it has been an internal face
            // of the original mesh and now it has become a face
            // on the parallel boundary
            // Use 0.5 weights
            addressing_[i].setSize(2);
            weights_[i].setSize(2, 0.5);

            addressing_[i][0] = own[ai];
            addressing_[i][1] = neighb[ai];
        }
        else
        {
            // This is a face that used to be on a cyclic boundary
            // but has now become a parallel patch face. I cannot
            // do the interpolation properly (I would need to look
            // up the different (face) list of data), so I will
            // just grab the value from the owner cell
            // HJ, 16/Mar/2001
            addressing_[i].setSize(1);
            weights_[i].setSize(1);

            addressing_[i][0] = own[ai];

            weights_[i][0] = 1.0;
        }
    }
}


fvFieldDecomposer::processorSurfacePatchFieldDecomposer::
processorSurfacePatchFieldDecomposer
(
    label sizeBeforeMapping,
    const unallocLabelList& addressingSlice,
    bool isScalarField
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    if (isScalarField)
    {
        // surfaceScalarFields represent fluxes. Thus, propogate
        // the correction of sign to account for flipped faces.
        forAll (addressing_, i)
        {
            addressing_[i].setSize(1);
            weights_[i].setSize(1);

            addressing_[i][0] = mag(addressingSlice[i]) - 1;
            weights_[i][0] = sign(addressingSlice[i]);
        }
    }
    else
    {
        // All surfaceFields other than of type scalar keep the sign unchanged
        forAll (addressing_, i)
        {
            addressing_[i].setSize(1);
            weights_[i].setSize(1);

            addressing_[i][0] = mag(addressingSlice[i]) - 1;
            weights_[i][0] = 1.0;
        }
    }
}


fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    boundaryAddressing_(boundaryAddressing),
    patchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<patchFieldDecomposer*>(nullptr)
    ),
    processorVolPatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorVolPatchFieldDecomposer*>(nullptr)
    ),
    processorSurfacePatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorSurfacePatchFieldDecomposer*>(nullptr)
    ),
    processorSurfaceScalarPatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorSurfacePatchFieldDecomposer*>(nullptr)
    )
{
    // HR 25.06.18: Weights cannot be used for mappings
    // HJ, 19/Oct/2018

    forAll (boundaryAddressing_, patchi)
    {
        if (boundaryAddressing_[patchi] >= 0)
        {
            patchFieldDecomposerPtrs_[patchi] = new patchFieldDecomposer
            (
                completeMesh_.boundary()[boundaryAddressing_[patchi]].size(),
                procMesh_.boundary()[patchi].patchSlice(faceAddressing_),
                completeMesh_.boundaryMesh()
                [
                    boundaryAddressing_[patchi]
                ].start()
            );
        }
        else
        {
            processorVolPatchFieldDecomposerPtrs_[patchi] =
                new processorVolPatchFieldDecomposer
                (
                    completeMesh_,
                    procMesh_.boundary()[patchi].patchSlice(faceAddressing_)
                );

            processorSurfacePatchFieldDecomposerPtrs_[patchi] =
                new processorSurfacePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    static_cast<const unallocLabelList&>
                    (
                        procMesh_.boundary()[patchi].patchSlice
                        (
                            faceAddressing_
                        )
                    ),
                    false
                );

            processorSurfaceScalarPatchFieldDecomposerPtrs_[patchi] =
                new processorSurfacePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    static_cast<const unallocLabelList&>
                    (
                        procMesh_.boundary()[patchi].patchSlice
                        (
                            faceAddressing_
                        )
                    ),
                    true
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fvFieldDecomposer::~fvFieldDecomposer()
{
    forAll (patchFieldDecomposerPtrs_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            delete patchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll (processorVolPatchFieldDecomposerPtrs_, patchi)
    {
        if (processorVolPatchFieldDecomposerPtrs_[patchi])
        {
            delete processorVolPatchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll (processorSurfacePatchFieldDecomposerPtrs_, patchi)
    {
        if (processorSurfacePatchFieldDecomposerPtrs_[patchi])
        {
            delete processorSurfacePatchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll (processorSurfaceScalarPatchFieldDecomposerPtrs_, patchi)
    {
        if (processorSurfaceScalarPatchFieldDecomposerPtrs_[patchi])
        {
            delete processorSurfaceScalarPatchFieldDecomposerPtrs_[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
