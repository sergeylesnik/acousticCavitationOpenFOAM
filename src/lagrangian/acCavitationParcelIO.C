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
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::acCavitationParcel::propHeader =
    KinematicParcel<acCavitationParcel>::propHeader
  + " dAv FBjPri";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Used for reading from file and for parallel transfer
Foam::acCavitationParcel::acCavitationParcel
(
    const Cloud<acCavitationParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    KinematicParcel<acCavitationParcel>(cloud, is, readFields),
    dAv_(0.0),
    VAv_(0.0),
    lastCell_(-1),
    nSmallDeltaStepFractions_(0),
    FBjPri_(vector::zero),
    PAc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            dAv_ = readScalar(is);
            is >> FBjPri_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&dAv_),
                sizeof(dAv_)
              + sizeof(FBjPri_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "acCavitationParcel::acCavitationParcel"
        "(const Cloud<acCavitationParcel>&, Istream&, bool)"
    );
}


void Foam::acCavitationParcel::readFields(Cloud<acCavitationParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    KinematicParcel<acCavitationParcel>::readFields(c);

    IOField<scalar> dAv(c.fieldIOobject("dAv", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dAv);

    IOField<vector> FBjPri(c.fieldIOobject("FBjPri", IOobject::MUST_READ));
    c.checkFieldIOobject(c, FBjPri);

    label i = 0;
    forAllIter(typename Cloud<acCavitationParcel>, c, iter)
    {
        acCavitationParcel& p = iter();

        p.dAv_ = dAv[i];
        p.FBjPri_ = FBjPri[i];
        i++;
    }
}


void Foam::acCavitationParcel::writeFields(const Cloud<acCavitationParcel>& c)
{
    KinematicParcel<acCavitationParcel>::writeFields(c);

    label np =  c.size();

    IOField<scalar> dAv(c.fieldIOobject("dAv", IOobject::NO_READ), np);
    IOField<vector> FBjPri(c.fieldIOobject("FBjPri", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<acCavitationParcel>, c, iter)
    {
        const acCavitationParcel& p = iter();

        dAv[i] = p.dAv();
        FBjPri[i] = p.FBjPri();
        i++;
    }

    dAv.write();
    FBjPri.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const acCavitationParcel& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const KinematicParcel<acCavitationParcel>&>(p)
            << token::SPACE << p.dAv()
            << token::SPACE << p.FBjPri();
    }
    else
    {
        os  << static_cast<const KinematicParcel<acCavitationParcel>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.dAv_),
            sizeof(p.dAv())
          + sizeof(p.FBjPri())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const acCavitationParcel&)"
    );

    return os;
}


// ************************************************************************* //
