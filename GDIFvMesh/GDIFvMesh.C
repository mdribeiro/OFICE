/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "GDIFvMesh.H"
#include "Time.H"
//#include "slidingInterface.H"
//#include "mapPolyMesh.H"
//#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H" 

#include "velocityLaplacianFvMotionSolver.H"
//#include "attachDetach.H"
#include "IFstream.H"
#include "interpolateXY.H" // for Liftprofile reading

#include "engineTime.H" // time convert to CAD and rad

#include "directionMixedFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GDIFvMesh, 0);

    addToRunTimeSelectionTable(dynamicFvMesh, GDIFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::GDIFvMesh::addZonesAndModifiers()
{
//#include   "addMeshModifierEngineGDIValve.H"
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::GDIFvMesh::GDIFvMesh(const IOobject& io)
:
    //topoChangerFvMesh(io),    
    dynamicFvMesh(io),
    //dynamicMeshDict
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    //GDI dictionary
    GDIdict_
    (
      IOdictionary
        (
	  IOobject
	  (
	      "GDIdict",
	      time().constant(),
	      *this,
	      IOobject::MUST_READ,
	      IOobject::NO_WRITE,
	      false
	  )
	)
    ),
    //Read engineGeometry just for Information at the very beginning   
    engineDict_
    (
        IOdictionary
        (
	  IOobject
	  (
            "engineGeometry",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
	  )
	)
    ),    
    //pointer for motionSolver
    msPtr_(motionSolver::New(*this)),
    
        ValveAxis_
    (
     (GDIdict_.lookup("ValveAxis"))
    ),
    
    // look for Engine Patch names
    piston_PatchName_
    (
        GDIdict_.lookup("piston_PatchName")
    ),
    //piston patch ID
    piston_PatchID_(-1),
   
   // Valve lift profile
    liftProfile1_
    (
        "theta",
        "lift",
        "Valve",
        IFstream
        (
           word(GDIdict_.lookup("liftProfileFile1"))
        )()
    ),
    
    liftProfile2_
    (
        "theta",
        "lift",
        "Valve",
        IFstream
        (
           word(GDIdict_.lookup("liftProfileFile2"))
        )()
    ),
    
    liftProfile3_
    (
        "theta",
        "lift",
        "Valve",
        IFstream
        (
           word(GDIdict_.lookup("liftProfileFile3"))
        )()
    ),
    liftProfile4_
    (
        "theta",
        "lift",
        "Valve",
        IFstream
        (
           word(GDIdict_.lookup("liftProfileFile4"))
        )()
    ),
    //Coordinate System for the Piston
   csPistonPtr_
    (
        coordinateSystem::New
        (
	   // piston_PatchName_,
            //"coordinateSystem",
            GDIdict_.subDict("Piston_CoordinateSystem")
        )
    ),
    //Coordinate System for the valve
    /*csValve1Ptr_
    (
        coordinateSystem::New
        (
            //"coordinateSystem",
            GDIdict_.subDict("Valve1_CoordinateSystem")
        )
    ),
    //Coordinate System for the valve
    csValve2Ptr_
    (
        coordinateSystem::New
        (
           // coordinateSystem,
            GDIdict_.subDict("Valve2_CoordinateSystem")
        )
    ),
    //Coordinate System for the valve
    csValve3Ptr_
    (
        coordinateSystem::New
        (
            //"coordinateSystem",
            GDIdict_.subDict("Valve3_CoordinateSystem")
        )
    ),
    //Coordinate System for the valve
    csValve4Ptr_
    (
        coordinateSystem::New
        (
           // "coordinateSystem",
            GDIdict_.subDict("Valve4_CoordinateSystem")
        )
    ),*/
    
     minLift_
    (
      readScalar(GDIdict_.lookup("minLift"))
    ),
    
    engineTime_(refCast<const engineTime>(time())),
    
     nValves_
    (
      readScalar(GDIdict_.lookup("nValves"))
    ),
     RPM_
    (
     (engineDict_.lookup("rpm"))
    ),
        CONROD_
    (
      (engineDict_.lookup("conRodLength"))
    ),
     BORE_
    (
     (engineDict_.lookup("bore"))
    ),
    STROKE_
    (
     (engineDict_.lookup("stroke"))
    ),
     CLEARANCE_
    (
     (engineDict_.lookup("clearance"))
    )

    
    

{
   Info <<"    ***ENGINE SPECIFICATION***" << nl;
   Info <<" Engine Speed     : " << RPM_.value()       << " rpm" << nl;
   Info <<" Engine Bore      : " << BORE_.value()      << " m" <<nl;
   Info <<" Engine Stroke    : " << STROKE_.value()    << " m" <<nl;
   Info <<" Engine Clearance : " << CLEARANCE_.value()    << " m" << " (not used within GDIFvMesh!) " <<nl; 
   Info <<"Number of Valves  :   "<<  nValves_ << nl;
   Info <<"Minimum Valve Lift: "<<  minLift_   << nl;
   Info <<"    ***READ PATCHES FOR MESH MOTION***" << nl;   piston_PatchID_ = boundaryMesh().findPatchID(piston_PatchName_);
   Info << " Piston ID: " <<  piston_PatchID_ << nl;
 /*  valveI =0;
    {
     // Valve  
     for(valveI=0;valveI< nValves_ ;valveI++)
     { 
     valveBottom_ =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_BottomPatchName"));
     valveStem_   =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_StemPatchName"  ));
     valveTop_    =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_TopPatchName"   ));
    
     Info << " Valve "<< valveI+1 <<" Bottom ID: " <<  valveBottom_ << nl;
     Info << " Valve "<< valveI+1 <<" Stem ID:   " <<  valveStem_ << nl;
     Info << " Valve "<< valveI+1 <<" Top ID:    " <<  valveTop_ << nl;
     }
    
    }*/
   //addZonesAndModifiers();
  // Pout << "Moving Profile:      Start CAD: " << min(liftProfile_.x()) << "    End CAD: " << max(liftProfile_.x()) <<nl;
   Info <<"    ***CONSTRUCTION FINISHED***" << nl;   piston_PatchID_ = boundaryMesh().findPatchID(piston_PatchName_);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GDIFvMesh::~GDIFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
  Foam::GDIFvMesh::csValve(int valve) const
 {
 
	switch (valve)
	      {
		case 1:
		  return csValve1Ptr_();
		case 2:
		  return csValve2Ptr_();
		case 3:
		  return csValve3Ptr_();
		case 4:
		  return csValve4Ptr_();
	      }
      }
*/
// interpolate Lift 

    Foam::scalar Foam::GDIFvMesh::lift(scalar theta ,int valve) const
    {
      
      switch (valve)
      {
	case 1:
	  return interpolateXY
	  (
	      theta,
	      liftProfile1_.x(),
	      liftProfile1_.y()
	  );
	case 2:
	 return interpolateXY
	  (
	      theta,
	      liftProfile2_.x(),
	      liftProfile2_.y()
	  );
	case 3:
	 return interpolateXY
	  (
	      theta,
	      liftProfile3_.x(),
	      liftProfile3_.y()
	  );
	case 4:
	 return interpolateXY
	  (
	      theta,
	      liftProfile4_.x(),
	      liftProfile4_.y()
	  );
      }
	  
    }
// current Lift 
    Foam::scalar Foam::GDIFvMesh::curLift(int valve) const
    {
	  return  lift(engineTime_.theta(), valve);
    }
    
// Valve speed
    Foam::scalar Foam::GDIFvMesh::curVel(int valve) const
    {
      if(isOpen(valve))
	  {
	  return (
	    ((curLift(valve)) - lift((engineTime_.theta() -  engineTime_.deltaTheta()), valve)
	    )/(engineTime_.deltaT().value() + VSMALL) //  or "time().deltaT().value()" it is the same!
	    );
	  }
	  else
	  {
	    return 0;
	  }
    }
    
//
bool Foam::GDIFvMesh::isOpen(int valve) const
{
    return lift(((engineTime_.theta())),valve) >= minLift_;
}
   
// Moving Mesh 
    bool Foam::GDIFvMesh::update()
    {
  
     // Info << "CAD =  " << engineTime_.theta() << nl;
      /*
      Info << "***  Interpolated Lift value:  CAD: " << engineTime_.theta() <<" lift: " << lift(engineTime_.theta()) << nl;
      Info << "     Valve Speed: "<< curVel() << " m/s" << nl;
    */  
     
      // Cast the msPtr_ to the demanded motion Solver. Here, we want to use the velocityLaplacianFvMotionSolver
      // but msPtr_() points only to fvMotionSolver
    
      velocityLaplacianFvMotionSolver& mSolver =
	    refCast<velocityLaplacianFvMotionSolver>(msPtr_());
	    
      #include  "setMotionBoundaryConditionGDIFvMesh.H"
	      
	// solve motion
	fvMesh::movePoints(msPtr_->newPoints());
    return 0; //opoChangeMap.valid();
    
}


// ************************************************************************* //
