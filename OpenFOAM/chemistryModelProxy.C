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

 ICE Revision: $Id: /local/openfoam/Libraries/alternateChemistryModels/OpenFOAM/chemistryModelProxy.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "chemistryModelProxy.H"
#include "addToRunTimeSelectionTable.H"
// #include "hMixtureThermo.H"

namespace Foam {

defineTypeNameAndDebug(chemistryModelProxy, 0);
addToRunTimeSelectionTable(alternateChemistryModel,chemistryModelProxy,steadyChemistry);
addToRunTimeSelectionTable(alternateChemistryModel,chemistryModelProxy,transientChemistry);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
chemistryModelProxy::chemistryModelProxy
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    alternateChemistryModel(thermo,rho),
    realChem_(thermo,rho)
{
//     if(!isA<hMixtureThermo<reactingMixture> >(thermo)){
//         FatalErrorIn("chemistryModelProxy::chemistryModelProxy")
//             << " thermo is required to be a hMixtureThermo<reactingMixture>"
//                 << endl << abort(FatalError);
//     }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

chemistryModelProxy::~chemistryModelProxy()
{}



// ************************************************************************* //
tmp<volScalarField> chemistryModelProxy::tc() const
{
    return realChem_.tc();
}

tmp<volScalarField> chemistryModelProxy::RR(const label i) const
{
    return realChem_.RR(i);
}

scalar chemistryModelProxy::solve(
    const scalar t0,
    const scalar deltaT
)
{
    return realChem_.solve(t0,deltaT);
}

void chemistryModelProxy::calcDQ(volScalarField &dQ)
{
    scalarField cp(dQ.size(), 0.0);
    const volScalarField &T=thermo_.T();

    forAll(realChem_.Y(), i)
    {
           volScalarField RRi = RR(i);

           forAll(T, celli)
           {
               scalar Ti = T[celli];
               cp[celli] += realChem_.Y()[i][celli]*realChem_.specieThermo()[i].Cp(Ti);
               scalar hi = realChem_.specieThermo()[i].h(Ti);
               scalar RR = RRi[celli];
               dQ[celli] -= hi*RR;
           }
           
    }
    
    forAll(dQ, celli)
    {
        dQ[celli] /= cp[celli];
    }     
}

} // namespace Foam
