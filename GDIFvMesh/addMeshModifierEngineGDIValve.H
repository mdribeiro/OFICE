{

    // Add a topology modifier
    Info << "Adding topology modifiers" << endl;

    PtrList<entry> entries( motionDict_.lookup("modifiers"));

   // Info << "Foam::name(valveI + 1)" << Foam::name(valveI + 1) << nl;
    //scalar a=1;
    topoChanger_.setSize(entries.size());
   
    forAll (entries, i)
    {
      Info << i << nl;
        topoChanger_.set
        (
            i,
            new attachDetach
            (
                entries[i].keyword(),
                i,
                topoChanger_,
	        //"detachFaceZoneV" + Foam::name(i + 1),
                word(entries[i].dict().lookup("attachFaceZoneV" + Foam::name(i + 1))),
	        word(entries[i].dict().lookup("valveDetachInCylV" + Foam::name(i + 1))),
	        word(entries[i].dict().lookup("valveDetachInPortV" + Foam::name(i+ 1))),
                scalarField(0),
                true                // Manual triggering
               // scalarField(entries[i].dict().lookup("triggerTimes"))
            )
        );
	Info << "Attach/Detach for Valve " << i << nl;
    }

    // Write mesh and modifiers
    
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}