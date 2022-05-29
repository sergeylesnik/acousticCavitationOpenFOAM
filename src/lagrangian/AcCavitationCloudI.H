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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::AcCavitationCloud<ParcelType>::constProps() const
{
    return constProps_;
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::AcCavitationCloud<ParcelType>::PAc() const
{
    return PAc_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::AcCavitationCloud<ParcelType>::argPAc() const
{
    return argPAc_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::AcCavitationCloud<ParcelType>::G() const
{
    return G_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::AcCavitationCloud<ParcelType>::argG() const
{
    return argG_;
}


template<class ParcelType>
inline bool Foam::AcCavitationCloud<ParcelType>::dampingViscous() const
{
    return dampingViscous_;
}


template<class ParcelType>
inline bool Foam::AcCavitationCloud<ParcelType>::dampingThermal() const
{
    return dampingThermal_;
}


template<class ParcelType>
inline bool Foam::AcCavitationCloud<ParcelType>::dampingRadiation() const
{
    return dampingRadiation_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::AcCavitationCloud<ParcelType>::PiViscous()
{
    return PiViscous_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::AcCavitationCloud<ParcelType>::PiThermal()
{
    return PiThermal_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::AcCavitationCloud<ParcelType>::PiRadiation()
{
    return PiRadiation_;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::AcCavitationCloud<ParcelType>::PiTotal
(
    Foam::volScalarField& PAc
)
{
    resetAcousticFields();
    
    tmp<volScalarField> tPiTotal
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "PiTotal",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
    );

    forAllConstIter(typename AcCavitationCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();
        const scalar REqu = p.REqu();
        const scalar PAcCell = PAc[cellI];

        if (dampingViscous_)
        {
            PiViscous_[cellI] += PiViTablePtr_()(REqu, PAcCell);
        }

        if (dampingThermal_)
        {
            PiThermal_[cellI] += PiThTablePtr_()(REqu, PAcCell);
        }

    }

    if (dampingViscous_)
    {
        tPiTotal().internalField() += PiViscous_;
    }

    if (dampingThermal_)
    {
        tPiTotal().internalField() += PiThermal_;
    }

    return tPiTotal;
}


template<class ParcelType>
inline Foam::acCavitationParticleForces&
Foam::AcCavitationCloud<ParcelType>::acCavForces()
{
    return acCavForces_;
}


template<class ParcelType>
inline Foam::label
Foam::AcCavitationCloud<ParcelType>::parcelsAddedConti() const
{
    return parcelsAddedConti_;
}


// ************************************************************************* //
