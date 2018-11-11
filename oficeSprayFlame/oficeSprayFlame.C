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
    sprayFoam

Description
    Transient PIMPLE solver for compressible, laminar or turbulent engine
flow swith spray parcels.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "engineTime.H"
#include "dynamicFvMesh.H"

#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiCombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "basicKinematicParcel.H"
#include "KinematicParcel.H"

#include "ignition.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    //#include "createTime.H"
    //#include "createMesh.H"
    #include "createEngineTime.H"
    #include "createDynamicFvMesh.H"

    #include "createControl.H"

   
    #include "createFields.H"
    #include "createFieldRefs.H"
 //   #include "createMRF.H"
    #include "createFvOptions.H"
    #include "createRhoUf.H"
    //#include "createClouds.H"
    //#include "createRadiationModel.H"
    #include "initContinuityErrs.H"
    #include "readEngineTimeControls.H"
    //#include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "startSummary.H"

    #include "readCombustionProperties.H"
  
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readEngineTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H" 

        runTime++;

        //Info<< "Time = " << runTime.timeName() << nl << endl;
	Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        //parcels.evolve();
	
	mesh.update(); 

	parcels.evolve();

        #include "rhoEqn.H"
	//#include "YEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            //#include "YEqn.H"
	    #include "ignite.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
		//#include "ignite.H"
		#include "YEqn.H"
	        #include "EEqn.H"
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        if (runTime.write())
        {
            combustion->Qdot()().write();
        }
        
        #include "logSummary.H"



	Info << "Max T = " << max(T).value() << " K      " << "Min T = " << min(T).value() << " K" << endl;
        Info << "Max p = " << max(p).value()/1.0e5 << " bar    " << "Min p = " << min(p).value()/1.0e5 << " bar" << endl;
        Info << "Mag(Max U) = " << mag(max(U)).value() << " m/s      " << "Mag(Min U) = " << mag(min(U)).value() << " m/s" << endl;
        Info << "Max C8H18 = " << max(composition.Y("C8H18")).value() << " --      " << "Min C8H18 = " << min(composition.Y("C8H18")).value() << " --" << endl;
	
	
        Info<< "ExecutionTiime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTiime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

