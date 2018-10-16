// Momentum predictor



// this is a hack to switch on or off the predictor or corrector time differencing shceme
// if predictorOrCorrector.value() < 0 , it is predictor step
// if predictorOrCorrector.value() > 0 , it is corrector step

{
    dimensionedScalar   predictorOrCorrector("predictorOrCorrector", 1);

    fvVectorMatrix UEqn
        (
            fvm::ddt(predictorOrCorrector,U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)-f
          
        );



    if (piso.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));     //Add Body force //TODO
    }

    // --- PISO loop
    while (piso.correct())
    {
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
                + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU);

        // Non-orthogonal pressure corrector loop
        while (piso.correctNonOrthogonal())
        {
            // Pressure corrector

            fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

            if (piso.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }                                                                                                                    

                #include "continuityErrs.H"

        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
    }

}
