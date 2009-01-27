/*---------------------------------------------------------------------------*\
This file written by Institute of Energy Process Enineering and Chemical
	Engineering TU Freiberg  http://www.iec.tu-freiberg.de
and ICE Stroemungsfoschungs GmbH http://www.ice-sf.at
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    alternateChemistryModel

Description

SourceFiles
    alternateChemistryModel.C

 ICE Revision: $Id: /local/openfoam/Libraries/alternateChemistryModels/alternateChemistryModel.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef alternateChemistryModel_H
#define alternateChemistryModel_H

#include "hCombustionThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class alternateChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class alternateChemistryModel
{

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const alternateChemistryModel&);

protected:
    
        const hCombustionThermo& thermo_;

        const volScalarField& rho_;

        tmp<volScalarField> characteristicTime();

public:

    TypeName("alternateChemistryModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        alternateChemistryModel,
        transientChemistry,
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        ),
        (thermo,rho)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        alternateChemistryModel,
        steadyChemistry,
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        ),
        (thermo,rho)
    );

    // Constructors

        //- Construct from components
        alternateChemistryModel
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        );


        static autoPtr<alternateChemistryModel> New
        (
            hCombustionThermo& thermo,
            const volScalarField& rho,
            bool steady=false
        );

    // Destructor

        virtual ~alternateChemistryModel();


    // Member Functions

        //- Return the thermodynamics model
        const hCombustionThermo& thermo() const
        {
            return thermo_;
        }

//         //- Access the fractions of the species
//         const PtrList<volScalarField>& Y() const
//         {
//             return Y_;
//         }

        //- Access the density field
        const volScalarField& rho() const 
        {
            return rho_;
        }

        //- The chemical time scale
        virtual tmp<volScalarField> tc() const=0;

	//- The local flow time scale
	tmp<volScalarField> tf();

        //- Return the chemical source terms
        virtual tmp<volScalarField> RR(const label i) const=0;

        //- Solve the reaction system for the given start time and time-step
        // and return the characteristic time
        virtual scalar solve
        (
            const scalar t0,
            const scalar deltaT
        ) = 0;

        //- Calculate dQ
        virtual void calcDQ(volScalarField &) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
