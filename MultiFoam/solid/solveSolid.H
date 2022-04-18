if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        Info << "\n \t solver solid hEqn !"<<endl;
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rho, h)
          - (
               thermo.isotropic()
             ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==            
            fvm::Su(HSource,h)
          + fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);

        Info << "\n \t solver solid FEqn !"<<endl;
        fvScalarMatrix FEqn
        (
            fvm::ddt(interF)
          - fvm::laplacian(Dsolid,interF)
        );
        FEqn.relax();
        fvOptions.constrain(FEqn);
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
