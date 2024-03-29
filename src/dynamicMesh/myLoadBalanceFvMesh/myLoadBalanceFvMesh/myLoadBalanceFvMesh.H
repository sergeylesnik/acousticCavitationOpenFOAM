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
    Foam::myLoadBalanceFvMesh

Description
    Load balancing fvMesh which is linked to a bug fixed decomposeReconstruct
    library

SourceFiles
    myLoadBalanceFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef myLoadBalanceFvMesh_H
#define myLoadBalanceFvMesh_H

#include "topoChangerFvMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class myLoadBalanceFvMesh Declaration
\*---------------------------------------------------------------------------*/

class myLoadBalanceFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dict_;

        //- Imbalance trigger 0 < imbalance < 1
        scalar imbalanceTrigger_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        myLoadBalanceFvMesh(const myLoadBalanceFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const myLoadBalanceFvMesh&);

        //- Check balance.  Return true if re-balancing is required
        bool checkLoadBalance(const scalarField& weights) const;


public:

    //- Runtime type information
    TypeName("myLoadBalanceFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit myLoadBalanceFvMesh(const IOobject& io);


    //- Destructor
    virtual ~myLoadBalanceFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
