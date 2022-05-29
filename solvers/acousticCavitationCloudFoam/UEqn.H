    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
        ==
        bubbleCloud.SU()/rho
    );

    if (pimple.momentumPredictor())
    {
        // Solve momentum predictor
        solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
    }
