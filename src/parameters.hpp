#ifndef pear_parameters_hpp
#define pear_parameters_hpp

#include <cmath>

namespace pear {

// SIMULATION PARAMETERS
const double Tcel = 25;
const double eta_u = 20.8 / 100;
const double eta_v = 0.04 / 100;

// OTHER AND DEPENDING PARAMETERS
const double T = Tcel + 273.15;
const double Tref = 293.15;
const double patm = 101300;
const double Rg = 8.314;

const double hu = 7e-7;
const double hv = 7.5e-7;

const double Cuamb = patm * eta_u / Rg / T;
const double Cvamb = patm * eta_v / Rg / T;

const double Vmu_ref = 2.39e-4;
const double Ea_vmu_ref = 80200;
const double Vmu = Vmu_ref * std::exp(Ea_vmu_ref / Rg * (1 / Tref - 1 / T));

const double Vmfv_ref = 1.61e-4;
const double Ea_vmfv_ref = 56700;
const double Vmfv = Vmfv_ref * std::exp(Ea_vmfv_ref / Rg * (1 / Tref - 1 / T));

const double Kmu = 0.4103;
const double Kmv = 27.2438;
const double Kmfu = 0.1149;
const double rq = 0.97;

const double Dur = 2.8e-10;
const double Duz = 1.1e-9;
const double Dvr = 2.32e-9;
const double Dvz = 6.97e-9;

} // namespace pear

#endif
