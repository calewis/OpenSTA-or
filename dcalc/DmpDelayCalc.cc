// OpenSTA, Static Timing Analyzer
// Copyright (c) 2026, Parallax Software, Inc.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
// 
// The origin of this software must not be misrepresented; you must not
// claim that you wrote the original software.
// 
// Altered source versions must be plainly marked as such, and must not be
// misrepresented as being the original software.
// 
// This notice may not be removed or altered from any source distribution.

#include "DmpDelayCalc.hh"

// ---- OpenROAD fork: two-pole load crossing solver (begin) ----
#include <algorithm>
#include <cmath>
// ---- OpenROAD fork: two-pole load crossing solver (end) ----
#include "DmpCeff.hh"
#include "GraphDelayCalc.hh"
#include "Liberty.hh"
#include "Network.hh"
#include "Parasitics.hh"
#include "PortDirection.hh"
#include "Sdc.hh"
#include "TableModel.hh"
#include "TimingArc.hh"

namespace sta {

// PiElmore parasitic delay calculator using Dartu/Menezes/Pileggi
// effective capacitance and elmore delay.
class DmpCeffElmoreDelayCalc : public DmpCeffDelayCalc
{
public:
  DmpCeffElmoreDelayCalc(StaState *sta);
  ArcDelayCalc *copy() override;
  std::string_view name() const override { return "dmp_ceff_elmore"; }
  ArcDcalcResult inputPortDelay(const Pin *port_pin,
                                float in_slew,
                                const RiseFall *rf,
                                const Parasitic *parasitic,
                                const LoadPinIndexMap &load_pin_index_map,
                                const Scene *scene,
                                const MinMax *min_max) override;

protected:
  void loadDelaySlew(const Pin *load_pin,
                     double drvr_slew,
                     const RiseFall *rf,
                     const LibertyLibrary *drvr_library,
                     const Parasitic *parasitic,
                     // Return values.
                     double &wire_delay,
                     double &load_slew) override;
};

ArcDelayCalc *
makeDmpCeffElmoreDelayCalc(StaState *sta)
{
  return new DmpCeffElmoreDelayCalc(sta);
}

DmpCeffElmoreDelayCalc::DmpCeffElmoreDelayCalc(StaState *sta) :
  DmpCeffDelayCalc(sta)
{
}

ArcDelayCalc *
DmpCeffElmoreDelayCalc::copy()
{
  return new DmpCeffElmoreDelayCalc(*this);
}

ArcDcalcResult
DmpCeffElmoreDelayCalc::inputPortDelay(const Pin *,
                                       float in_slew,
                                       const RiseFall *rf,
                                       const Parasitic *parasitic,
                                       const LoadPinIndexMap &load_pin_index_map,
                                       const Scene *scene,
                                       const MinMax *min_max)
{
  Parasitics *parasitics = scene->parasitics(min_max);
  ArcDcalcResult dcalc_result(load_pin_index_map.size());
  LibertyLibrary *drvr_library = network_->defaultLibertyLibrary();
  for (auto [load_pin, load_idx] : load_pin_index_map) {
    double wire_delay = 0.0;
    double load_slew = in_slew;
    bool elmore_exists = false;
    float elmore = 0.0;
    if (parasitic)
      parasitics->findElmore(parasitic, load_pin, elmore, elmore_exists);
    if (elmore_exists)
      // Input port with no external driver.
      dspfWireDelaySlew(load_pin, rf, in_slew, elmore, wire_delay, load_slew);
    thresholdAdjust(load_pin, drvr_library, rf, wire_delay, load_slew);
    dcalc_result.setWireDelay(load_idx, wire_delay);
    dcalc_result.setLoadSlew(load_idx, load_slew);
  }
  return dcalc_result;
}

void
DmpCeffElmoreDelayCalc::loadDelaySlew(const Pin *load_pin,
                                      double drvr_slew,
                                      const RiseFall *rf,
                                      const LibertyLibrary *drvr_library,
                                      const Parasitic *parasitic,
                                      // Return values.
                                      double &wire_delay,
                                      double &load_slew)
{
  wire_delay = 0.0;
  load_slew = drvr_slew;
  bool elmore_exists = false;
  float elmore = 0.0;
  if (parasitic)
    parasitics_->findElmore(parasitic, load_pin, elmore, elmore_exists);
  if (elmore_exists) {
    if (auto r = loadDelaySlewElmore(load_pin, elmore)) {
      wire_delay = r->first;
      load_slew = r->second;
    }
  }
  thresholdAdjust(load_pin, drvr_library, rf, wire_delay, load_slew);
}

////////////////////////////////////////////////////////////////

// PiPoleResidue parasitic delay calculator using Dartu/Menezes/Pileggi
// effective capacitance and two poles/residues.
class DmpCeffTwoPoleDelayCalc : public DmpCeffDelayCalc
{
public:
  DmpCeffTwoPoleDelayCalc(StaState *sta);
  ArcDelayCalc *copy() override;
  std::string_view name() const override { return "dmp_ceff_two_pole"; }
  Parasitic *findParasitic(const Pin *drvr_pin,
                           const RiseFall *rf,
                           const Scene *scene,
                           const MinMax *min_max) override;
  Parasitic *reduceParasitic(const Parasitic *parasitic_network,
                             const Pin *drvr_pin,
                             const RiseFall *rf,
                             const Scene *scene,
                             const MinMax *min_max) override;
  ArcDcalcResult inputPortDelay(const Pin *port_pin,
                                float in_slew,
                                const RiseFall *rf,
                                const Parasitic *parasitic,
                                const LoadPinIndexMap &load_pin_index_map,
                                const Scene *scene,
                                const MinMax *min_max) override;
  ArcDcalcResult gateDelay(const Pin *drvr_pin,
                           const TimingArc *arc,
                           const Slew &in_slew,
                           float load_cap,
                           const Parasitic *parasitic,
                           const LoadPinIndexMap &load_pin_index_map,
                           const Scene *scene,
                           const MinMax *min_max) override;

protected:
  void loadDelaySlew(const Pin *load_pin,
                     double drvr_slew,
                     const RiseFall *rf,
                     const LibertyLibrary *drvr_library,
                     const Parasitic *parasitic,
                     // Return values.
                     double &wire_delay,
                     double &load_slew) override;

  using ArcDelayCalc::reduceParasitic;

private:
  void loadDelay(double drvr_slew,
                 Parasitic *pole_residue,
                 double p1,
                 double k1,
                 double &wire_delay,
                 double &load_slew);
  // ---- OpenROAD fork: two-pole load crossing solver (begin) ----
  static double twoPoleLoadCrossing(double vth,
                                    double p1,
                                    double p2,
                                    double k1,
                                    double k2,
                                    double tt,
                                    double y_tt);
  // ---- OpenROAD fork: two-pole load crossing solver (end) ----

  bool parasitic_is_pole_residue_{false};
  float vth_{0.0F};
  float vl_{0.0F};
  float vh_{0.0F};
  float slew_derate_{0.0F};
};

ArcDelayCalc *
makeDmpCeffTwoPoleDelayCalc(StaState *sta)
{
  return new DmpCeffTwoPoleDelayCalc(sta);
}

DmpCeffTwoPoleDelayCalc::DmpCeffTwoPoleDelayCalc(StaState *sta) :
  DmpCeffDelayCalc(sta)
{
}

ArcDelayCalc *
DmpCeffTwoPoleDelayCalc::copy()
{
  return new DmpCeffTwoPoleDelayCalc(*this);
}

Parasitic *
DmpCeffTwoPoleDelayCalc::findParasitic(const Pin *drvr_pin,
                                       const RiseFall *rf,
                                       const Scene *scene,
                                       const MinMax *min_max)
{
  Parasitic *parasitic = nullptr;
  const Sdc *sdc = scene->sdc();
  Parasitics *parasitics = scene->parasitics(min_max);
  if (parasitics == nullptr
      // set_load net has precedence over parasitics.
      || network_->direction(drvr_pin)->isInternal())
    return nullptr;
  // Prefer PiPoleResidue.
  parasitic = parasitics->findPiPoleResidue(drvr_pin, rf, min_max);
  if (parasitic)
    return parasitic;
  parasitic = parasitics->findPiElmore(drvr_pin, rf, min_max);
  if (parasitic)
    return parasitic;
  Parasitic *parasitic_network =
    parasitics->findParasiticNetwork(drvr_pin);
  if (parasitic_network) {
    parasitic = reduceParasitic(parasitic_network, drvr_pin, rf,
                                scene, min_max);
    if (parasitic)
      return parasitic;
  }
  Wireload *wireload = sdc->wireload(min_max);
  if (wireload) {
    float pin_cap, wire_cap, fanout;
    bool has_wire_cap;
    graph_delay_calc_->netCaps(drvr_pin, rf, scene, min_max, pin_cap, wire_cap,
                               fanout, has_wire_cap);
    parasitic = parasitics->estimatePiElmore(drvr_pin, rf, wireload,
                                             fanout, pin_cap,
                                             scene, min_max);
  }
  return parasitic;
}

// Reduce to PiPoleResidue, the format loadDelaySlew consumes; the
// inherited PiElmore reduction yields zero wire delays.
Parasitic *
DmpCeffTwoPoleDelayCalc::reduceParasitic(const Parasitic *parasitic_network,
                                         const Pin *drvr_pin,
                                         const RiseFall *rf,
                                         const Scene *scene,
                                         const MinMax *min_max)
{
  Parasitics *parasitics = scene->parasitics(min_max);
  return parasitics->reduceToPiPoleResidue2(parasitic_network, drvr_pin, rf,
                                            scene, min_max);
}

ArcDcalcResult
DmpCeffTwoPoleDelayCalc::inputPortDelay(const Pin *,
                                        float in_slew,
                                        const RiseFall *rf,
                                        const Parasitic *parasitic,
                                        const LoadPinIndexMap &load_pin_index_map,
                                        const Scene *scene,
                const MinMax *min_max)
{
  const Parasitics *parasitics = scene->parasitics(min_max);
  ArcDcalcResult dcalc_result(load_pin_index_map.size());
  double wire_delay = 0.0;
  double load_slew = in_slew;
  LibertyLibrary *drvr_library = network_->defaultLibertyLibrary();
  for (const auto [load_pin, load_idx] : load_pin_index_map) {
    if (parasitics->isPiPoleResidue(parasitic)) {
      const Parasitic *pole_residue = parasitics->findPoleResidue(parasitic, load_pin);
      if (pole_residue) {
        size_t pole_count = parasitics->poleResidueCount(pole_residue);
        if (pole_count >= 1) {
          ComplexFloat pole1, residue1;
          // Find the 1st (elmore) pole.
          parasitics->poleResidue(pole_residue, 0, pole1, residue1);
          if (pole1.imag() == 0.0
              && residue1.imag() == 0.0) {
            float p1 = pole1.real();
            float elmore = 1.0F / p1;
            dspfWireDelaySlew(load_pin, rf, in_slew, elmore, wire_delay, load_slew);
            thresholdAdjust(load_pin, drvr_library, rf, wire_delay, load_slew);
          }
        }
      }
    }
    dcalc_result.setWireDelay(load_idx, wire_delay);
    dcalc_result.setLoadSlew(load_idx, load_slew);
  }
  return dcalc_result;
}

ArcDcalcResult
DmpCeffTwoPoleDelayCalc::gateDelay(const Pin *drvr_pin,
                                   const TimingArc *arc,
                                   const Slew &in_slew,
                                   float load_cap,
                                   const Parasitic *parasitic,
                                   const LoadPinIndexMap &load_pin_index_map,
                                   const Scene *scene,
                                   const MinMax *min_max)
{
  parasitics_ = scene->parasitics(min_max);
  const LibertyLibrary *drvr_library = arc->to()->libertyLibrary();
  const RiseFall *rf = arc->toEdge()->asRiseFall();
  vth_ = drvr_library->outputThreshold(rf);
  vl_ = drvr_library->slewLowerThreshold(rf);
  vh_ = drvr_library->slewUpperThreshold(rf);
  slew_derate_ = drvr_library->slewDerateFromLibrary();
  return DmpCeffDelayCalc::gateDelay(drvr_pin, arc, in_slew, load_cap, parasitic,
                                     load_pin_index_map, scene, min_max) ;
}

void
DmpCeffTwoPoleDelayCalc::loadDelaySlew(const Pin *load_pin,
                                       double drvr_slew,
                                       const RiseFall *rf,
                                       const LibertyLibrary *drvr_library,
                                       const Parasitic *parasitic,
                                       // Return values.
                                       double &wire_delay,
                                       double &load_slew)
{
  parasitic_is_pole_residue_ = parasitics_->isPiPoleResidue(parasitic);
  // ---- OpenROAD fork: PiElmore load delay fallback (begin) ----
  wire_delay = 0.0;
  load_slew = drvr_slew;
  if (!parasitic_is_pole_residue_) {
    if (parasitics_->isPiElmore(parasitic)) {
      float elmore = 0.0F;
      bool exists = false;
      parasitics_->findElmore(parasitic, load_pin, elmore, exists);
      if (exists) {
        if (auto elmore_delay_slew = loadDelaySlewElmore(load_pin, elmore)) {
          wire_delay = elmore_delay_slew->first;
          load_slew = elmore_delay_slew->second;
        }
      }
    }
    thresholdAdjust(load_pin, drvr_library, rf, wire_delay, load_slew);
    return;
  }
  Parasitic *pole_residue = parasitics_->findPoleResidue(parasitic, load_pin);
  // ---- OpenROAD fork: PiElmore load delay fallback (end) ----
  if (pole_residue) {
    size_t pole_count = parasitics_->poleResidueCount(pole_residue);
    if (pole_count >= 1) {
      ComplexFloat pole1, residue1;
      // Find the 1st (elmore) pole.
      parasitics_->poleResidue(pole_residue, 0, pole1, residue1);
      if (pole1.imag() == 0.0
          && residue1.imag() == 0.0) {
        float p1 = pole1.real();
        float k1 = residue1.real();
        if (pole_count >= 2)
          loadDelay(drvr_slew, pole_residue, p1, k1, wire_delay, load_slew);
        // ---- OpenROAD fork: single-pole ramp load crossing (begin) ----
        else if (p1 > 0.0F) {
          if (!delayZero(drvr_slew, this)) {
            double tt = drvr_slew * slew_derate_ / (vh_ - vl_);
            double p1_d = p1;
            double y_tt = (tt + std::expm1(-p1_d * tt) / p1_d) / tt;
            wire_delay = twoPoleLoadCrossing(vth_, p1_d, 0.0, p1_d, 0.0, tt, y_tt)
                - tt * vth_;
            double tl = twoPoleLoadCrossing(vl_, p1_d, 0.0, p1_d, 0.0, tt, y_tt);
            double th = twoPoleLoadCrossing(vh_, p1_d, 0.0, p1_d, 0.0, tt, y_tt);
            load_slew = std::max(0.0, (th - tl) / slew_derate_);
          }
          else {
            wire_delay = -std::log(1.0 - vth_) / p1;
            load_slew = std::log((1.0 - vl_) / (1.0 - vh_)) / (p1 * slew_derate_);
          }
        }
        // ---- OpenROAD fork: single-pole ramp load crossing (end) ----
      }
    }
  }
  thresholdAdjust(load_pin, drvr_library, rf, wire_delay, load_slew);
}

void
DmpCeffTwoPoleDelayCalc::loadDelay(double drvr_slew,
                                   Parasitic *pole_residue,
                                   double p1,
                                   double k1,
                                   // Return values.
                                   double &wire_delay,
                                   double &load_slew)
{
  ComplexFloat pole2, residue2;
  parasitics_->poleResidue(pole_residue, 1, pole2, residue2);
  if (!delayZero(drvr_slew, this)
      && pole2.imag() == 0.0
      && residue2.imag() == 0.0) {
    double p2 = pole2.real();
    double k2 = residue2.real();
    // ---- OpenROAD fork: two-pole load crossing solver (begin) ----
    if (!(p1 > 0.0) || !(p2 > 0.0) || !std::isfinite(k1) || !std::isfinite(k2)) {
      wire_delay = p1;
      load_slew = p1;
      return;
    }
    double k1_p1_2 = k1 / (p1 * p1);
    double k2_p2_2 = k2 / (p2 * p2);
    // Convert tt to 0:1 range.
    double tt = drvr_slew * slew_derate_ / (vh_ - vl_);
    double y_tt =
        (tt + k1_p1_2 * std::expm1(-p1 * tt) + k2_p2_2 * std::expm1(-p2 * tt)) / tt;
    wire_delay = twoPoleLoadCrossing(vth_, p1, p2, k1, k2, tt, y_tt) - tt * vth_;

    double tl = twoPoleLoadCrossing(vl_, p1, p2, k1, k2, tt, y_tt);
    double th = twoPoleLoadCrossing(vh_, p1, p2, k1, k2, tt, y_tt);
    load_slew = std::max(0.0, (th - tl) / slew_derate_);
    // ---- OpenROAD fork: two-pole load crossing solver (end) ----
  }
}

// ---- OpenROAD fork: two-pole load crossing solver (begin) ----
double
DmpCeffTwoPoleDelayCalc::twoPoleLoadCrossing(double vth,
                                             double p1,
                                             double p2,
                                             double k1,
                                             double k2,
                                             double tt,
                                             double y_tt)
{
  const double k1_p1_2 = k1 / (p1 * p1);
  const double k2_p2_2 = (k2 != 0.0) ? k2 / (p2 * p2) : 0.0;
  const double c1_post = -k1_p1_2 * std::expm1(-p1 * tt);
  const double c2_post = (k2 != 0.0) ? -k2_p2_2 * std::expm1(-p2 * tt) : 0.0;

  auto eval = [&](double t, double &f, double &df) {
    if (t <= tt) {
      double y =
          (t + k1_p1_2 * std::expm1(-p1 * t) + k2_p2_2 * std::expm1(-p2 * t)) / tt;
      double yp = (1.0 - (k1 / p1) * std::exp(-p1 * t)
                   - ((k2 != 0.0) ? (k2 / p2) * std::exp(-p2 * t) : 0.0))
          / tt;
      f = y - vth;
      df = yp;
    }
    else {
      double u = t - tt;
      double e1 = std::exp(-p1 * u);
      double e2 = (k2 != 0.0) ? std::exp(-p2 * u) : 0.0;
      double y = 1.0 - (c1_post * e1 + c2_post * e2) / tt;
      double yp = (p1 * c1_post * e1 + p2 * c2_post * e2) / tt;
      f = y - vth;
      df = yp;
    }
  };

  double lo, hi, t;
  if (y_tt < vth) {
    lo = tt;
    double p_min = (k2 != 0.0) ? std::min(p1, p2) : p1;
    double span = std::max(tt, 1.0 / p_min);
    hi = tt + span;
    double f_hi = 0.0;
    double df_hi = 0.0;
    for (int i = 0; i < 64; ++i) {
      eval(hi, f_hi, df_hi);
      if (f_hi >= 0.0)
        break;
      lo = hi;
      span *= 2.0;
      hi = tt + span;
    }
    double denom = (1.0 - vth) * p1 * p1 * tt;
    double log_arg =
        (denom > 0.0 && k1 > 0.0) ? (k1 * (-std::expm1(-p1 * tt))) / denom : 0.0;
    t = (log_arg > 0.0) ? tt + std::log(log_arg) / p1 : 0.5 * (lo + hi);
    if (!(t > lo && t < hi))
      t = 0.5 * (lo + hi);
  }
  else {
    lo = 0.0;
    hi = tt;
    t = (y_tt > 0.0) ? vth * tt / y_tt : 0.5 * tt;
    if (!(t > lo && t < hi))
      t = 0.5 * (lo + hi);
  }

  for (int iter = 0; iter < 64; ++iter) {
    double f = 0.0;
    double df = 0.0;
    eval(t, f, df);
    if (std::abs(f) <= 1e-13 || (hi - lo) <= 1e-13 * std::max(1.0, hi))
      return t;
    if (f < 0.0)
      lo = t;
    else
      hi = t;
    double t_next = (df > 0.0) ? (t - f / df) : 0.5 * (lo + hi);
    if (!(t_next > lo && t_next < hi))
      t_next = 0.5 * (lo + hi);
    if (std::abs(t_next - t) <= 1e-13 * std::max(1.0, t_next))
      return t_next;
    t = t_next;
  }
  return t;
}
// ---- OpenROAD fork: two-pole load crossing solver (end) ----

} // namespace sta
