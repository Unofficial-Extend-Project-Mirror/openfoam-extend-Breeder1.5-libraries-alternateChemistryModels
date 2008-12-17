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

 ICE Revision: $Id: /local/openfoam/Libraries/alternateChemistryModels/OpenFOAM/localTimeChemistryModelProxy.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "localTimeChemistryModelProxy.H"
#include "addToRunTimeSelectionTable.H"

#include "chemistrySolver.H"

namespace Foam {

defineTypeNameAndDebug(localTimeChemistryModelProxy, 0);
addToRunTimeSelectionTable(alternateChemistryModel,localTimeChemistryModelProxy , steadyChemistry);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
localTimeChemistryModelProxy::localTimeChemistryModelProxy
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    chemistryModelProxy(thermo,rho)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

localTimeChemistryModelProxy::~localTimeChemistryModelProxy()
{}



// ************************************************************************* //

scalar localTimeChemistryModelProxy::solve(
    const scalar t0,
    const scalar defaultDeltaT
)
{
    for(label i=0; i<realChem_.Ns(); i++)
    {
        realChem_.RR_[i].setSize(rho_.size());
    }

    if (!realChem_.chemistry_)
    {
        return GREAT;
    }

    scalar deltaTMin = GREAT;
    const scalarField localTime=characteristicTime()().internalField();

    for(label celli=0; celli<rho_.size(); celli++)
    {
        scalar deltaT=min(defaultDeltaT,localTime[celli]);

        for(label i=0; i<realChem_.Ns(); i++)
        {
            realChem_.RR_[i][celli] = 0.0;
        }

        scalar rhoi = rho_[celli];
        scalar Ti = thermo_.T()[celli];
        scalar hi = thermo_.h()[celli];
        scalar pi = thermo_.p()[celli];

        scalarField c(realChem_.Ns_);
        scalarField c0(realChem_.Ns_);
        scalarField dc(realChem_.Ns_, 0.0);

        for(label i=0; i<realChem_.Ns_; i++)
        {
            c[i] = rhoi*realChem_.Y_[i][celli]/realChem_.specieThermo_[i].W();
        }
        c0 = c;

        scalar t = t0;
        scalar tauC = realChem_.deltaTChem_[celli];
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        // calculate the chemical source terms
        scalar cTot = 0.0;

        while(timeLeft > SMALL)
        {
            tauC = realChem_.solver().solve(c, Ti, pi, t, dt);
            t += dt;

            // update the temperature
            cTot = sum(c);
            chemistryModel::reactionThermo mixture(0.0*realChem_.specieThermo_[0]);
            for(label i=0; i<realChem_.Ns_; i++)
            {
                mixture += (c[i]/cTot)*realChem_.specieThermo_[i];
            }        
            Ti = mixture.TH(hi, Ti);

            timeLeft -= dt;
            realChem_.deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;
        scalar WTot = 0.0;
        for(label i=0; i<realChem_.Ns_; i++)
        {
            WTot += c[i]*realChem_.specieThermo_[i].W();
        }        
        WTot /= cTot;

        for(label i=0; i<realChem_.Ns_; i++)
        {
            realChem_.RR_[i][celli] = dc[i]*realChem_.specieThermo_[i].W()/deltaT;
        }
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*defaultDeltaT);

    return deltaTMin;
}

} // namespace Foam
