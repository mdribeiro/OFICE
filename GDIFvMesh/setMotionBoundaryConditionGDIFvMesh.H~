// set Boundary velocity for the solid 
Info << "Setting Boundary Conditions for motion solver" << endl;
valveI =0;
{
    pointVectorField& pointMotionU =  mSolver.pointMotionU();
    
    // Valve  
     for(valveI=0;valveI< nValves_ ;valveI++)
   
     //forAll (topoChanger_, valveI)
     {
     Info << "ValveAxis_" << ValveAxis_[valveI] << endl;  
     valveBottom_ =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_BottomPatchName"));
     valveStem_   =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_StemPatchName"  ));
     valveTop_    =  boundaryMesh().findPatchID(GDIdict_.lookup("valve" + Foam::name(valveI+1) + "_TopPatchName"   ));  
              
     pointMotionU.boundaryField()[valveBottom_] ==  scalar(-1) * curVel(valveI+1) *  ValveAxis_[valveI];//csValve(valveI+1).axis();	
     pointMotionU.boundaryField()[valveStem_]   ==  scalar(-1) * curVel(valveI+1) *  ValveAxis_[valveI];//csValve(valveI+1).axis();	
     pointMotionU.boundaryField()[valveTop_]    ==  scalar(-1) * curVel(valveI+1) *  ValveAxis_[valveI];//csValve(valveI+1).axis();
     
     Info << "Valve "<< valveI+1 <<": "<< "Velocity: "<<  curVel(valveI+1)<<"  Lift:" << lift(engineTime_.theta(),valveI+1) << nl;

       
       // Info << "current Valve Lift for Valve "<< valveI <<": "<< lift(engineTime_.theta(),valveI)<< nl;
     };
     Info << "Piston Position: "<< engineTime_.pistonPosition().value() << " m" << "  Piston Speed: "<< engineTime_.pistonSpeed().value() << " m/s"<< nl;
  //  Info << "b.c for valves set" << nl;
  //  Info <<"finished piston position"<< csPistonPtr_.Rtr_.e3()<<endl;
     
     vector axis ;
     axis =  vector ( 0,0,1);
// Piston    
     Info << "Piston Position: "<<piston_PatchID_<< endl;
    pointMotionU.boundaryField()[piston_PatchID_] == axis* engineTime_.pistonSpeed().value();
// Cell motion B.C. for liner 

       volVectorField& cellMotionU =  mSolver.cellMotionU();
// find liner patch      
           label patchID = boundaryMesh().findPatchID("liner");
	    if(patchID == -1)
	    {
		Info << "patch not found" << endl;
		return 0;
	    }
	    else
	    {
	    directionMixedFvPatchVectorField& UcellMotion =
		    refCast<directionMixedFvPatchVectorField>
		    (
		    cellMotionU.boundaryField()[patchID]
		    );
// liner 
	    //- patch face normals      
	    vectorField n = boundary()[patchID].nf();
	   // Pout << "face normals of liner" << n << nl;
	    
	    //- set value fraction so as to have fixed value in patch normal direction
	    UcellMotion.valueFraction() = symm(n*n);

	    //- set fixed value
	    UcellMotion.refValue() = vector (0,0, 0); 
	    //- set gradient
            UcellMotion.refGrad() = vector (0, 0, 0);
          
	    }
      //  Info << "b.c for liner set" << nl;
    }
Info << "Boundary Conditions for motion solver are set" << endl;


