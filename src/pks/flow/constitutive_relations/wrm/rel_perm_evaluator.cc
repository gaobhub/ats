/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "rel_perm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    min_val_(0.) {

  ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms),
    min_val_(0.) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrms_(other.wrms_),
    sat_key_(other.sat_key_),
    min_val_(other.min_val_) {}


Teuchos::RCP<FieldEvaluator>
RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


void RelPermEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<string>("rel perm key", "relative_permeability");
  }
  setLinePrefix(my_key_+std::string(" evaluator"));

  // my dependencies are just saturation.
  sat_key_ = plist_.get<string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum rel perm cutoff", 1.e-30);
}


void RelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  if (!wrms_->first->initialized()) wrms_->first->Initialize(result->mesh());

  const Epetra_MultiVector& sat = *S->GetFieldData(sat_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& res_v = *result->ViewComponent("cell",false);

  // Evaluate the evaluator to calculate sat.
  int ncells = res_v.MyLength();
  for (int c=0; c!=ncells; ++c) {
    int index = (*wrms_->first)[c];
    double pc = wrms_->second[index]->capillaryPressure(sat[0][c]);
    res_v[0][c] = std::max(wrms_->second[index]->k_relative(pc), min_val_);
  }
}


void RelPermEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  // not implemented, not yet needed --etc
  ASSERT(0);
}



} //namespace
} //namespace
} //namespace
