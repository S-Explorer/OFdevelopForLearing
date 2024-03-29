if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        Info << "\n \t solver solid FEqn !"<<endl;
        fvScalarMatrix FEqn
        (
            fvm::ddt(interF)
          - fvm::laplacian(Dsolid,interF)
        );
        FEqn.relax();
        fvOptions.constrain(FEqn);
        FEqn.solve();
        
		Info << "\n \t solver solid hEqn !"<<endl;
        fvScalarMatrix hEqn
        (
            fvm::ddt(rho, h)
            // - fvm::ddt(rho(1-interF*cp_g/cp), h)
          - fvm::laplacian(kappa/cp, h, "laplacian(alpha,h)")
          - fvc::ddt(rho * cp/ cp_g, interF)
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);

        FEqn.solve();
    }

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
