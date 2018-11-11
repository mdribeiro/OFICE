/*if (ign.ignited())
{
      #include "EEqn.H"
}*/



  forAll(ign.sites(), i)
    {
        const ignitionSite& ignSite = ign.sites()[i];

        if (ignSite.igniting())
        {
            forAll(ignSite.cells(), icelli)
            {
                label ignCell = ignSite.cells()[icelli];
              //  Info<< "Igniting cell " << ignCell;
		   //const volScalarField& T = thermo.T();  
		   // 2411
		   T[ignSite.cells()[icelli]] = Tig.value(); // Temperature of adiabatic flame for gasoline

		   Info<< " Hi ignition! Tig = " << Tig.value() << endl;

                /*  Info<< " state :"
                    << ' ' << c[ignCell]
                    << ' ' << Xi[ignCell]
                    << ' ' << Su[ignCell]
                    << ' ' << mgb[ignCell]
                    << endl; */
		    //c[ignSite.cells()[icelli]] = 1;
		   // Info<< " c after  :" <<  c[ignCell] << endl;
		    //Tb[ignSite.cells()[icelli]] = 2500;
		    //w[ignSite.cells()[icelli]] =  ignSite.strength() ; //8000;

            }
        }
    }
