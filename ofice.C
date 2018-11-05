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
    engineFoam

Description
    Solver for cold-flow in internal combustion engines.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
//#include "engineMesh.H"
#include "dynamicFvMesh.H"

#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "OFstream.H"
#include "fvIOoptionList.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createEngineTime.H"
    //#include "createEngineMesh.H"
    #include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "createRhoUf.H"
    #include "initContinuityErrs.H"
    #include "readEngineTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "startSummary.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readEngineTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"   
      
        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg"
            << endl;

        //mesh.move();
	     mesh.update(); 

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
	    //#include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "EEqn.H"
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        
	Info << "Max T = " << max(T).value() << " K      " << "Min T = " << min(T).value() << " K" << endl;
	Info << "Max rho = " << max(rho).value() << " kg/m3      " << "Min rho = " << min(rho).value() << " kg/m3" << endl;
        Info << "Max p = " << max(p).value()/1.0e5 << " bar    " << "Min p = " << min(p).value()/1.0e5 << " bar" << endl;
	Info << "Mag(Max U) = " << mag(max(U)).value() << " m/s      " << "Mag(Min U) = " << mag(min(U)).value() << " m/s" << endl;        
      
/*
       volScalarField C(sqrt(thermo.Cp()/thermo.Cv()*(1.0/psi)));


	volScalarField Mach
	(
	    IOobject
	    (
		"Mach",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mag(U)/mag(C)
	);	
	
	Info << "Mag(Max C) = " << mag(max(C)).value() << " m/s      " << "Mag(Min C) = " << mag(min(C)).value() << " m/s" << endl; 
	Info << "Max Mach = " << max(Mach).value() << " -      " << "Min Mach = " << min(Mach).value() << " -" << endl; 	
	
	*/

        runTime.write();

        #include "logSummary.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
