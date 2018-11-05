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

Class
    Foam::GDIFvMesh

Description


SourceFiles
    GDIFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef GDIFvMesh_H
#define GDIFvMesh_H

//#include "topoChangerFvMesh.H"
#include "dynamicFvMesh.H"
#include "motionSolver.H"

#include "graph.H"
#include "word.H"

#include "coordinateSystem.H"
#include "engineTime.H" // time convert to CAD and rad
#include "directionMixedFvPatchFields.H"

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "cylindricalCS.H"

#include "vectorList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class GDIFvMesh Declaration
\*---------------------------------------------------------------------------*/

class GDIFvMesh
:
    public dynamicFvMesh //topoChangerFvMesh //, public engineTime
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;
	
	//- GDI dictionary
        dictionary GDIdict_;

	//- engineGeometry dictionary
        dictionary engineDict_;
	
        //- Motion solver
        autoPtr<motionSolver> msPtr_;
		
        //patch names to move in GDIdict
	word piston_PatchName_;
	//- piston patch ID
        label piston_PatchID_;	

	//- Lift profile
        graph liftProfile1_;
		
	//- Lift profile
        graph liftProfile2_;
	
	//- Lift profile
        graph liftProfile3_;
	
	//- Lift profile
        graph liftProfile4_;
	
	//- Coordinate system Piston
        autoPtr<coordinateSystem> csPistonPtr_;
	
	//- Coordinate system
        autoPtr<coordinateSystem> csValve1Ptr_;
	
	//- Coordinate system
        autoPtr<coordinateSystem> csValve2Ptr_;
	
	//- Coordinate system
        autoPtr<coordinateSystem> csValve3Ptr_;
	
	//- Coordinate system
        autoPtr<coordinateSystem> csValve4Ptr_;

	//- Minimum valve lift.  On this lift the valve is considered closed
         scalar minLift_;
	
	 // numbering for attach/Detach
	int valveI;
	// counter
	//int valve;

	label valveBottom_;
	label valveStem_;
	label valveTop_;
	 
	scalar nValves_;
	
	//Engine Specification Info
	dimensionedScalar RPM_;
	dimensionedScalar CONROD_;
	dimensionedScalar BORE_;
	dimensionedScalar STROKE_;
	dimensionedScalar CLEARANCE_;
	
	//Axis for Valves
	
	vectorList ValveAxis_;
	
    // Private Member Functions

        //- Disallow default bitwise copy construct
        GDIFvMesh(const GDIFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const GDIFvMesh&);

 //- Add attach/detach zones and modifiers
        void addZonesAndModifiers();
	
	
protected:
	//engine Database
	const engineTime& engineTime_;
  
public:

    //- Runtime type information
    TypeName("GDIFvMesh");


    // Constructors

        //- Construct from database
        explicit GDIFvMesh(const IOobject& io);


    // Destructor

        virtual ~GDIFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

	//- Return coordinate system
        const coordinateSystem& csPiston() const
        {
            return csPistonPtr_();
        }
        
        const coordinateSystem& csValve(int valve) const
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
               
        const engineTime& engTime() const
        {
            return engineTime_;
        }
	
    // solid position and velocity

	  //- Return valve lift given in time space
            scalar lift(scalar theta, int valve)  const; 
	    
	 // Return current lift position
	    scalar curLift(int valve) const;
	    
	 // Return current lift velocity
	    scalar curVel(int valve) const ;
	  
	 //- Is the valve open?
         bool isOpen(int valve) const;
	    
	    
	   
	    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
