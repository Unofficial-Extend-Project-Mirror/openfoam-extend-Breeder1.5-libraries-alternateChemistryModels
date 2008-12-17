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
    localTimeChemistryModelProxy

Description

SourceFiles
    localTimeChemistryModelProxy.C

 ICE Revision: $Id: /local/openfoam/Libraries/alternateChemistryModels/OpenFOAM/localTimeChemistryModelProxy.H 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#ifndef localTimeChemistryModelProxy_H
#define localTimeChemistryModelProxy_H

// don't try this at home kids!!!
#define private public

#include "chemistryModelProxy.H"

// Setting it back to normal
#undef private

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localTimeChemistryModelProxy Declaration
\*---------------------------------------------------------------------------*/

class localTimeChemistryModelProxy
:
    public chemistryModelProxy
{

        //- Disallow default bitwise assignment
        void operator=(const localTimeChemistryModelProxy&);

public:

    TypeName("localTimeChemistryModelProxy");

    // Constructors

        //- Construct from components
        localTimeChemistryModelProxy
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        );


    // Destructor

        virtual ~localTimeChemistryModelProxy();

        //- Solve the reaction system for the given start time and time-step
        // and return the characteristic time
        virtual scalar solve
        (
            const scalar t0,
            const scalar deltaT
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
