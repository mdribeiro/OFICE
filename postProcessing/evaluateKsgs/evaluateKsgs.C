/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    evaluateKsgs

Description
    Utility for the calculation of the ratio of sub-grid k and total k (sub-grid + resolved)
    turbRatio = 0 means DNS
    turbRatio = 1 means RANS

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"

#include "fvCFD.H"
#include "engineTime.H"
#include "dynamicFvMesh.H"

#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"

#include "argList.H"
#include "OFstream.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createEngineTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    Foam::timeSelector::addOptions();
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info << nl << "Time = " << runTime.timeName() << nl << endl;

	Info << "Max T = " << max(T).value() << " K      " << "Min T = " << min(T).value() << " K" << endl;
        Info << "Max p = " << max(p).value()/1.0e5 << " bar    " << "Min p = " << min(p).value()/1.0e5 << " bar" << endl;
        Info << "Mag(Max U) = " << mag(max(U)).value() << " m/s      " << "Mag(Min U) = " << mag(min(U)).value() << " m/s" << endl;

	Info << "Get k_sgs!" << endl;
	ksgs_ = turbulence->k();

	scalar count = 0;
	scalar U_avg = 0;
	forAll(U,celli)
        {
                U_avg = U_avg + mag(U[celli]);
		count = count + 1;
        }
	U_avg = U_avg/count;

	forAll(mesh.C(),celli)
	{
		Kres[celli] = 0.5*( mag(U[celli]) - U_avg  )*( mag(U[celli]) - U_avg  );
	}

	forAll(mesh.C(),celli)
        {
                turbRate[celli] = ksgs_[celli]/( Kres[celli] + ksgs_[celli]  );
        }


	count = 0;
        scalar turbRate_avg = 0;
        forAll(U,celli)
        {
                turbRate_avg = turbRate_avg + turbRate[celli];
                count = count + 1;
        }
        turbRate_avg = turbRate_avg/count;

	ksgs_.write();
	Kres.write();
	turbRate.write();

	Info << "turbRate avg = " << turbRate_avg << endl;

        Info<< "ExecutionTiime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTiime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
