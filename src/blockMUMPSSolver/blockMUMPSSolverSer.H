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
    blockMUMPSSolverSer

Description
    An implementation of the abstract blockMUMPSSolver for the solution of
    linear systems with the direct solver MUMPS in a serial run.

Author
    Sergey Lesnik, ITM Clausthal

SourceFiles
    blockMUMPSSolverSer.C

\*---------------------------------------------------------------------------*/

#ifndef blockMUMPSSolverSer_H
#define blockMUMPSSolverSer_H

#include "blockMUMPSSolver.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blockMUMPSSolverSer Declaration
\*---------------------------------------------------------------------------*/

class blockMUMPSSolverSer
:
    public blockMUMPSSolver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        blockMUMPSSolverSer(const blockMUMPSSolverSer&);

        //- Disallow default bitwise assignment
        void operator=(const blockMUMPSSolverSer&);


public:

    //- Runtime type information
    TypeName("blockMUMPSSolverSer")

    // Constructors

        //- Construct from matrix components and mesh
        blockMUMPSSolverSer
        (
            fvBlockMatrix<vector2>& matrix,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~blockMUMPSSolverSer()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        void solveCore();

        //- Get global index for MUMPS of the diagonal entries
        inline label diagInd(label i)
        {
            return i*2;
        }

        //- Set global indices for MUMPS of the owner and neighbour for a face
        void setOwnerNeighbInd
        (
            label &ownerI,
            label &neighbourI,
            const label &faceI
        );

        //- Dump current complete RHS
        void dumpCompleteRhs(fileName filePrefix);

        //- Dump according to MUMPS dictionary from case
        void dumpAccToMUMPSDict();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
