/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#include "Mesh_Algorithms.hh"
#include "pk_helpers.hh"

namespace Amanzi {

bool
aliasVector(State& S, const Key& key, const Tag& target, const Tag& alias)
{
  if (S.HasEvaluator(key, target) && !S.HasEvaluator(key, alias)) {
    S.SetEvaluator(key, alias, S.GetEvaluatorPtr(key, target));
    S.GetRecordSetW(key).AliasRecord(target, alias);
    return true;
  }
  return false;
}


// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector.
// -----------------------------------------------------------------------------
void
applyDirichletBCs(const Operators::BCs& bcs, CompositeVector& u)
{
  if (u.HasComponent("face")) {
    Epetra_MultiVector& u_f = *u.ViewComponent("face",false);
    for (unsigned int f=0; f!=u_f.MyLength(); ++f) {
      if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_f[0][f] = bcs.bc_value()[f];
      }
    }
  }

  if (u.HasComponent("boundary_face")) {
    Epetra_MultiVector& u_bf = *u.ViewComponent("boundary_face", false);
    const Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
    const Epetra_Map& vandalay_map = u.Mesh()->exterior_face_map(false);
    const Epetra_Map& face_map = u.Mesh()->face_map(false);

    for (int bf=0; bf!=u_bf.MyLength(); ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = bcs.bc_value()[f];
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Given a vector and a face ID, get the value at that location.
//
// Looks in the following order:
//  -- face component
//  -- boundary Dirichlet data
//  -- boundary_face value
//  -- internal cell
// -----------------------------------------------------------------------------
double
getFaceOnBoundaryValue(AmanziMesh::Entity_ID f, const CompositeVector& u, const Operators::BCs& bcs)
{
  if (u.HasComponent("face")) {
    return (*u.ViewComponent("face",false))[0][f];
  } else if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
    return bcs.bc_value()[f];
  // } else if (u.HasComponent("boundary_face")) {
  //   AmanziMesh::Entity_ID bf = getFaceOnBoundaryBoundaryFace(*u.Mesh(), f);
  //   return (*u.ViewComponent("boundary_face",false))[0][bf];
  } else {
    auto c = getFaceOnBoundaryInternalCell(*u.Mesh(),f);
    return (*u.ViewComponent("cell",false))[0][c];
  }
  return -1;
}


// -----------------------------------------------------------------------------
// Get the directional int for a face that is on the boundary.
// -----------------------------------------------------------------------------
int
getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f)
{
  AmanziMesh::Entity_ID_List cells;
  mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  AMANZI_ASSERT(cells.size() == 1);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh.cell_get_faces_and_dirs(cells[0], &faces, &dirs);
  return dirs[std::find(faces.begin(), faces.end(), f) - faces.begin()];
}


// -----------------------------------------------------------------------------
// Get a primary variable evaluator for a key at tag
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
RequireEvaluatorPrimary(const Key& key, const Tag& tag, State& S)
{
  // first check, is there one already
  if (S.HasEvaluator(key, tag)) {
    // if so, make sure it is primary
    Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
    Teuchos::RCP<EvaluatorPrimaryCV> eval_pv =
      Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
    if (eval_pv == Teuchos::null) {
      Errors::Message msg;
      msg << "Expected primary variable evaluator for "
          << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }
    return eval_pv;
  }

  // if not, create one, only at this tag, not to be shared across tags.  By
  // this, we mean we don't stick the "type" = "primary" back into the
  // evaluator list -- this allows "copy evaluators" e.g. "water content at the
  // old tag" to differ from the standard evalulator, e.g. "water content at
  // the new tag" which is likely a secondary variable evaluator.
  Teuchos::ParameterList plist(key);
  plist.set("evaluator type", "primary variable");
  plist.set("tag", tag.get());
  auto eval_pv = Teuchos::rcp(new EvaluatorPrimaryCV(plist));
  S.SetEvaluator(key, tag, eval_pv);
  return eval_pv;
}


// -----------------------------------------------------------------------------
// Marks a primary evaluator as changed.
// -----------------------------------------------------------------------------
void
ChangedEvaluatorPrimary(const Key& key, const Tag& tag, State& S)
{
  Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
  Teuchos::RCP<EvaluatorPrimaryCV> eval_pv =
    Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
  if (eval_pv == Teuchos::null) {
    Errors::Message msg;
    msg << "Expected primary variable evaluator for "
        << key << " @ " << tag.get();
    Exceptions::amanzi_throw(msg);
  }
  eval_pv->SetChanged();
}


} // namespace Amanzi
