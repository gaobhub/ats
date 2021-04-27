/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell with subgrid microtopography.

/*!

  Uses the subgrid equation from Jan et al WRR 2018 for volumetric or
  effective ponded depth to determine the area of water, then heuristically
  places snow on top of that surface.

Requires the following dependencies:

* `"snow depth key`" ``[string]`` **DOMAIN-depth**
         The name of the snow depth.

* `"snow transition height [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the snow gradient.
* `"water transition height [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the water gradient.
* `"minimum fractional area [-]`" ``[double]`` **1.e-5**
         Mimimum area fraction allowed, less than this is rebalanced as zero.

Ordering of the area fractions calculated are: [bare ground, water, snow].

NOTE: this evaluator simplifies the situation by assuming constant
density.  This make it so that ice and water see the same geometry per
unit pressure, which isn't quite true thanks to density differences.
However, we hypothesize that these differences, on the surface (unlike in
the subsurface) really don't matter much. --etc

*/

#include "area_fractions_threechannel_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// Constructor from ParameterList
AreaFractionsThreeChannelEvaluator::AreaFractionsThreeChannelEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  snow_subgrid_transition_ = plist_.get<double>("snow transition height [m]", 0.02);
  water_subgrid_transition_ = plist_.get<double>("water transition height [m]", 0.02);
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ <= 0.) {
    Errors::Message message("AreaFractionsThreeChannelEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  domain_ = Keys::getDomain(my_key_);
  if (domain_ == "surface") {
    domain_snow_ = "snow";
  } else if (Keys::starts_with(domain_, "surface_")) {
    domain_snow_ = std::string("snow_") + domain_.substr(8,domain_.size());
  }

  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
  ponded_depth_key_ = Keys::readKey(plist_, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(ponded_depth_key_);
}


void
AreaFractionsThreeChannelEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto& res = *result->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);
  const auto& pd = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell",false);

  for (int c=0; c!=res.MyLength(); ++c) {
    // calculate area of land
    if (sd[0][c] >= snow_subgrid_transition_) {
      res[1][c] = 0.;
      res[2][c] = 1.;
    } else {
      if (sd[0][c] <= 0.) {
        res[2][c] = 0.;
      } else {
        res[2][c] = sd[0][c] / snow_subgrid_transition_;
      }

      // snow preferentially covers water, as both go to low lying areas
      if (pd[0][c] >= water_subgrid_transition_) {
        res[1][c] = 1 - res[2][c];
      } else if (pd[0][c] <= 0.) {
        res[1][c] = 0.;
      } else {
        double water_covered = pd[0][c] / water_subgrid_transition_;
        if (res[2][c] > water_covered) {
          res[1][c] = 0;
        } else {
          res[1][c] = water_covered - res[2][c];
        }
      }
    }
    res[0][c] = 1 - res[1][c] - res[2][c];

    // if any area is less than eps, give to others
    // if any area fraction is less than eps, give it to the others
    if (res[0][c] > 0 && res[0][c] < min_area_) {
      if (res[1][c] < min_area_) {
        res[2][c] = 1.;
        res[1][c] = 0.;
        res[0][c] = 0.;
      } else {
        res[1][c] += res[0][c] * res[1][c] / (res[1][c] + res[2][c]);
        res[2][c] += res[0][c] * res[2][c] / (res[1][c] + res[2][c]);
        res[0][c] = 0.;
      }
    } else if (res[1][c] > 0 && res[1][c] < min_area_) {
      if (res[2][c] < min_area_) {
        res[0][c] = 1.;
        res[1][c] = 0.;
        res[2][c] = 0.;
      } else {
        res[0][c] += res[1][c] * res[0][c] / (res[0][c] + res[2][c]);
        res[2][c] += res[1][c] * res[2][c] / (res[0][c] + res[2][c]);
        res[1][c] = 0.;
      }
    } else if (res[2][c] > 0 && res[2][c] < min_area_) {
      res[0][c] += res[2][c] * res[0][c] / (res[0][c] + res[1][c]);
      res[1][c] += res[2][c] * res[1][c] / (res[0][c] + res[1][c]);
      res[2][c] = 0.;
    }

    AMANZI_ASSERT(std::abs(res[0][c] + res[1][c] + res[2][c] - 1.0) < 1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[0][c] && res[0][c] <= 1.+1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[1][c] && res[1][c] <= 1.+1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[2][c] && res[1][c] <= 1.+1.e-10);

    res[0][c] = std::min(std::max(0.,res[0][c]), 1.);
    res[1][c] = std::min(std::max(0.,res[1][c]), 1.);
    res[2][c] = std::min(std::max(0.,res[2][c]), 1.);
  }
}

void
AreaFractionsThreeChannelEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // see if we can find a master fac
  auto my_fac = S->RequireField(my_key_, my_key_);
  my_fac->SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Check plist for vis or checkpointing control.
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  for (auto dep_key : dependencies_) {
    auto fac = S->RequireField(dep_key);
    if (Keys::getDomain(dep_key) == domain_snow_) {
      fac->SetMesh(S->GetMesh(domain_snow_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else {
      fac->SetMesh(S->GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
    }

    // Recurse into the tree to propagate info to leaves.
    S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace

