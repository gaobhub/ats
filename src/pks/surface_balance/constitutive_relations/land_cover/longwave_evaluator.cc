/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates incoming longwave radiation from rel humidity and air temperature.

/*!

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec

    * `"minimum relative humidity [-]`" ``[double]`` **0.1** Sets a minimum rel humidity, RH=0 breaks the model.

    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"relative humidity key`" ``[string]`` **DOMAIN-relative_humidity**

*/

#include "Key.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "longwave_evaluator.hh"


namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

LongwaveEvaluator::LongwaveEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  air_temp_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{air_temp_key_, tag});
  rel_hum_key_ = Keys::readKey(plist, domain, "relative humidity", "relative_humidity");
  dependencies_.insert(KeyTag{rel_hum_key_, tag});

  min_rel_hum_ = plist.get<double>("minimum relative humidity [-]", 0.1);
  scale_ = plist.get<double>("scaling factor [-]", 1.0);
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
LongwaveEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& air_temp = *S.Get<CompositeVector>(air_temp_key_, tag).ViewComponent("cell", false);
  const auto& rel_hum = *S.Get<CompositeVector>(rel_hum_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    res[0][c] = scale_ * Relations::IncomingLongwaveRadiation(air_temp[0][c], std::max(min_rel_hum_, rel_hum[0][c]));
  }
}

} //namespace
} //namespace
} //namespace

