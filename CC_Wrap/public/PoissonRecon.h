#ifndef POISSON_RECON_INCLUDED
#define POISSON_RECON_INCLUDED
template<typename Real>
struct PoissonReconInput{
    Real position[3];
};
template<typename Real>
struct PoissonReconOutput{
    Real position[3];
};
template<typename Real>
PoissonReconOutput<Real> ExecutePoissonRecon(PoissonReconInput<Real> input);
#endif