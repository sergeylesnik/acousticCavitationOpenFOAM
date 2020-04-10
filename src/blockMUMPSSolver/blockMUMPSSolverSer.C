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

#include "blockMUMPSSolverSer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(blockMUMPSSolverSer, 0)
    addToRunTimeSelectionTable
    (
        blockMUMPSSolver,
        blockMUMPSSolverSer,
        runType
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMUMPSSolverSer::blockMUMPSSolverSer
(
    fvBlockMatrix<vector2> &matrix,
    const fvMesh &mesh
)
:
    blockMUMPSSolver(matrix, mesh)
{
    // Diagonal is asSqaure 2x2 -> 2* entries in solution vector
    // and order of the matrix.
    mumps_.n = 2*nCells_;

    // Matrix is non-symmetric for MUMPS since the diagonal
    // coupling terms have different signs -> 4*NfaceB;
    // diagonal is 2x2 -> 4* entries in the matrix A
    mumps_.nz = 4*nInternalFaces_ + 4*nCells_;

    irn_[myid_].resize(mumps_.nz);
    jcn_[myid_].resize(mumps_.nz);
    amv_[myid_].resize(mumps_.nz);
    rhs_[myid_].resize(mumps_.n);

    assembleDiag();
    assembleOffDiag();

    mumps_.irn = irn_[myid_].data();
    mumps_.jcn = jcn_[myid_].data();
    mumps_.a = amv_[myid_].data();

    analyzeAndFactorizeMUMPS();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::blockMUMPSSolverSer::solveCore()
{
    mumps_.rhs = rhs_[myid_].data();
    solveMUMPS();
}

void Foam::blockMUMPSSolverSer::setOwnerNeighbInd
(
    label &ownerI,
    label &neighbourI,
    const label &faceI
)
{
    // lowerAddr are the indices of the owners of the face
    ownerI = 2 * matrix_.lduAddr().lowerAddr()[faceI];

    // upperAddr are the indices of the neighbours of the owner cells
    // which share the same face
    neighbourI = 2 * matrix_.lduAddr().upperAddr()[faceI];
}

void Foam::blockMUMPSSolverSer::dumpCompleteRhs(fileName filePrefix)
{
    OFstream osRhs
    (
        filePrefix + "RHS_Time" + mesh_.time().timeName()
    );

    osRhs << "procNo localI globRow rhs"<< endl;

    forAll(rhs_, procI)
    {
        forAll(rhs_[procI], i)
        {
            osRhs.precision(12);
            osRhs << procI << " " << i
                  << " " << i+1  // Fortran style
                  << " " << rhs_[procI][i] << endl;
        }
    }
}

void Foam::blockMUMPSSolverSer::dumpAccToMUMPSDict()
{
    if (dumpCompleteLinSys_)
    {
        dumpMasterLinearSystem(dumpFilePrefix_);
    }
}
