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

 ICE Revision: $Id: /local/openfoam/Libraries/alternateChemistryModels/alternateChemistryModel.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

#include "alternateChemistryModel.H"

#include "chemistryModelProxy.H"
#include "hMixtureThermo.H"

namespace Foam {

defineTypeNameAndDebug(alternateChemistryModel, 0);

defineRunTimeSelectionTable(alternateChemistryModel, transientChemistry);
defineRunTimeSelectionTable(alternateChemistryModel, steadyChemistry);

// Construct from components
alternateChemistryModel::alternateChemistryModel
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    thermo_(thermo),
    rho_(rho)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

alternateChemistryModel::~alternateChemistryModel()
{}

autoPtr<alternateChemistryModel> alternateChemistryModel::New
(
    hCombustionThermo& thermo,
    const volScalarField& rho,
    bool steady
)
{
    Info << "Choosing the chemistry engine" << endl;

    word chemistryEngineName;

    {
        IOdictionary chemistryPropertiesDict
        (
            IOobject
            (
                "chemistryProperties",
                rho.time().constant(),
                rho.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        chemistryPropertiesDict.lookup("chemistryEngine")
            >> chemistryEngineName;
    }

    Info<< "Selecting chemistry engine " << chemistryEngineName << endl;

    if(steady) {
        steadyChemistryConstructorTable::iterator cstrIter =
            steadyChemistryConstructorTablePtr_->find(chemistryEngineName);

        if (cstrIter == steadyChemistryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "alternateChemistryModel::New"
                )   << "Unknown chemistryEngine type " << chemistryEngineName
                    << endl << endl
                    << "Valid steady chemistry engines are :" << endl
                    << steadyChemistryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }
	Info<< "with STEADY chemistry model" << endl;

        return autoPtr<alternateChemistryModel>(cstrIter()(thermo,rho));
    } else {
        transientChemistryConstructorTable::iterator cstrIter =
            transientChemistryConstructorTablePtr_->find(chemistryEngineName);

        if (cstrIter == transientChemistryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "alternateChemistryModel::New"
                )   << "Unknown chemistryEngine type " << chemistryEngineName
                    << endl << endl
                    << "Valid transient chemistry engines are :" << endl
                    << transientChemistryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }
	Info<< "with UNSTEADY chemistry model" << endl;

        return autoPtr<alternateChemistryModel>(cstrIter()(thermo,rho));
    }
}

tmp<volScalarField> alternateChemistryModel::characteristicTime()
{
    autoPtr<volVectorField> Utmp;
    if(!rho_.mesh().foundObject<volVectorField>("U")) {
        WarningIn("alternateChemistryModel::characteristicTime()")
            << "No velocity field found .... force reading" << endl;

        Utmp.reset(
            new volVectorField 
            (
                IOobject
                (
                    "U",
                    rho_.mesh().time().timeName(),
                    rho_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                rho_.mesh()
            )
        );

    }

    const volVectorField &U=rho_.mesh().lookupObject<volVectorField>("U");
    
    volScalarField charLen
        (
            IOobject
            (
                "charLen",
                rho_.mesh().time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho_.mesh(),
            dimensionedScalar("L",dimensionSet(0,1,0,0,0,0,0),0)
        );
    
        // #define LENGTH_BY_VOLUME
#ifdef LENGTH_BY_VOLUME
    charLen.internalField()=pow(rho_.mesh().V(),1./3);
#else
    const pointField &pts=rho_.mesh().points();
    const cellList &cells=rho_.mesh().cells();
    const faceList &faces=rho_.mesh().faces();
    
    forAll(cells,cellI) {
        if(mag(U[cellI])>SMALL) {
            vector dir=U[cellI]/mag(U[cellI]);
            scalar maxX=-GREAT;
            scalar minX= GREAT;
            const cell &c=cells[cellI];
            forAll(c,fI) {
                const face &f=faces[c[fI]];
                forAll(f,vI) {
                    scalar x= dir & pts[f[vI]];
                    if(x>maxX) {
                        maxX=x;
                    }
                    if(x<minX) {
                        minX=x;
                    }
                }
            }
            charLen[cellI]=maxX-minX;
        } else {
            charLen[cellI]=GREAT;
        }
    }
#endif

    tmp<volScalarField> charTime(
        new volScalarField(
            IOobject
            (
                "charTime",
                rho_.mesh().time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho_.mesh(),
            dimensionedScalar("T",dimensionSet(0,0,1,0,0,0,0),0)
        )
    );

    forAll(charTime(),cellI) {
        if(mag(U[cellI])>SMALL) {
            charTime()[cellI]=charLen[cellI]/mag(U[cellI]);
        } else {
            charTime()[cellI]=GREAT;
        }
    }

    return charTime;
}

tmp<volScalarField> alternateChemistryModel::tf()
{
	const scalarField localTime=characteristicTime()().internalField();

	tmp<volScalarField> tsource(
	    new volScalarField(
		 IOobject
		 (
		      "tf",
                      rho_.time().timeName(),
		      rho_.db(),
		      IOobject::NO_READ,
		      IOobject::NO_WRITE
		 ),
		 rho_.mesh(),
		 dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0),
		 zeroGradientFvPatchScalarField::typeName
		 )
	    );

	tsource().internalField() = localTime;
	tsource().correctBoundaryConditions();

	return tsource;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

}

// ************************************************************************* //
