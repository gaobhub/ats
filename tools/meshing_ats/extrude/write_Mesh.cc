#include <set>
#include <vector>
#include <algorithm>
#include "exodusII.h"

#include "dbc.hh"

#include "write_Mesh.hh"


namespace Amanzi {
namespace AmanziGeometry {

void
writeExodus(const Mesh3D& m, const std::string& filename) {
  // create the exodus file
  int CPU_word_size = sizeof(float);
  int IO_Word_size = 8;
  int fid = ex_create(filename.c_str(), EX_NOCLOBBER, &CPU_word_size, &IO_Word_size);
  if (fid < 0) {
    std::cerr << "Cowardly not clobbering: \"" << filename << "\" already exists." << std::endl;
    return;
  }
  
  // make the blocks by set
  std::set<int> set_ids(m.block_ids.begin(), m.block_ids.end());
  std::vector<std::vector<int> > blocks;
  std::vector<int> blocks_ncells;
  std::vector<int> blocks_id;
  std::vector<std::vector<int> > block_face_counts;
  
  for (auto sid : set_ids) {
    std::vector<int> block;
    std::vector<int> face_counts;
    // count things
    int ncells = 0;

    for (int i=0; i!=m.cell2face.size(); ++i) {
      if (m.block_ids[i] == sid) {
        ncells++;
        block.insert(block.end(), m.cell2face[i].begin(), m.cell2face[i].end());
        face_counts.push_back(m.cell2face[i].size());
      }
    }
    

    blocks.emplace_back(block);
    blocks_ncells.push_back(ncells);
    blocks_id.push_back(sid);
    block_face_counts.emplace_back(face_counts);
  }


  // create the params
  ex_init_params params;
  sprintf(params.title, "my_mesh");
  params.num_dim = 3;
  params.num_nodes = m.coords.size();
  params.num_edge = 0;
  params.num_edge_blk = 0;
  params.num_face = m.face2node.size();
  params.num_face_blk = 1;
  params.num_elem = m.cell2face.size();
  params.num_elem_blk = set_ids.size();
  params.num_node_maps = 0;
  params.num_edge_maps = 0;
  params.num_face_maps = 0;
  params.num_elem_maps = 0;
  params.num_side_sets = m.face_sets.size();
  params.num_elem_sets = 0;

  int ierr = ex_put_init_ext(fid, &params);
  ASSERT(!ierr);
  
  // make the coordinate arrays, set the coordinates
  // NOTE: exodus seems to only deal with floats!
  std::vector<std::vector<float> > coords(3);
  for (int i=0; i!=3; ++i) {
    coords[i].resize(m.coords.size());
  }

  for (int n=0; n!=coords[0].size(); ++n) {
    coords[0][n] = m.coords[n][0];
    coords[1][n] = m.coords[n][1];
    coords[2][n] = m.coords[n][2];
  }

  char* coord_names[3];
  char a[10]="xcoord";
  char b[10]="ycoord";
  char c[10]="zcoord";
  coord_names[0]=a;
  coord_names[1]=b;
  coord_names[2]=c;

  ierr |= ex_put_coord_names(fid, coord_names);
  ASSERT(!ierr);
  ierr |= ex_put_coord(fid, &coords[0][0], &coords[1][0], &coords[2][0]);
  ASSERT(!ierr);

  // put in the face block
  std::vector<int> facenodes;
  std::vector<int> facenodes_counts;
  for (auto& nodes : m.face2node) {
    facenodes_counts.push_back(nodes.size());
    facenodes.insert(facenodes.end(), nodes.begin(), nodes.end());
  }
  ierr |= ex_put_block(fid, EX_FACE_BLOCK, 1, "NSIDED",
                       m.face2node.size(), facenodes.size(), 0,0,0);
  ASSERT(!ierr);


  ierr |= ex_put_entity_count_per_polyhedra(fid, EX_FACE_BLOCK, 1,
          &facenodes_counts[0]);
  ASSERT(!ierr);

  for (auto& e : facenodes) e++;
  ierr |= ex_put_conn(fid, EX_FACE_BLOCK, 1, NULL, NULL, &facenodes[0]);
  ASSERT(!ierr);
  

  // put in the element blocks
  for (int lcvb=0; lcvb!=blocks.size(); ++lcvb) {
    ierr |= ex_put_block(fid, EX_ELEM_BLOCK, blocks_id[lcvb], "NFACED",
                         blocks_ncells[lcvb], 0, 0, blocks[lcvb].size(),0);
    ASSERT(!ierr);

    ierr |= ex_put_entity_count_per_polyhedra(fid, EX_ELEM_BLOCK, blocks_id[lcvb],
            &block_face_counts[lcvb][0]);
    ASSERT(!ierr);

    for (auto&e : blocks[lcvb]) e++;
    ierr |= ex_put_conn(fid, EX_ELEM_BLOCK, blocks_id[lcvb], NULL, NULL, &blocks[lcvb][0]);
    ASSERT(!ierr);
  }

  // add the side sets
  for (int lcvs=0; lcvs!=m.face_sets.size(); ++lcvs) {
    auto& s = m.face_sets[lcvs];
    auto elems_copy(s.first);
    auto faces_copy(s.second);
    for (auto& e : elems_copy) e++;
    for (auto& e : faces_copy) e++;
    ierr |= ex_put_side_set_param(fid, m.face_sets_id[lcvs], elems_copy.size(), 0);
    ierr |= ex_put_side_set(fid, m.face_sets_id[lcvs], &elems_copy[0], &faces_copy[0]);
  }


  ierr |= ex_close(fid);

}

}
}
