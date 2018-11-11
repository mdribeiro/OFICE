// Solve the Momentum equation

//MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
 ==
    rho()*g
  + parcels.SU(U)
  + fvOptions(rho, U)
);

//UEqn().relax();

//fvOptions.constrain(UEqn());

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
