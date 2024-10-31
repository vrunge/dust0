#include "2D_B_Indices2.h"
#include <cmath>
using namespace Rcpp;

// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //


Indices_2D2::Indices_2D2(){}

Indices_2D2::~Indices_2D2(){}

void Indices_2D2::reset() { current = list.begin(); }
void Indices_2D2::next() { ++current; }
bool Indices_2D2::check() { return current != list.end(); }

void Indices_2D2::set_init_size(const unsigned int& size) { list.reserve(size); }
void Indices_2D2::add(const unsigned int& value){list.push_back(value);}

std::vector<unsigned int> Indices_2D2::get_list() { return list; }
void Indices_2D2::remove_last() { list.pop_back(); }



////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// ////  DeterministicIndices_2D  ////
///// ///// ///// ///// ///// ///// /////

DeterministicIndices_2D::DeterministicIndices_2D() : Indices_2D2() {}

// full reset for pruning step
void DeterministicIndices_2D::reset_prune()
{
  // reset iterators
  if (list.size() > 1)
  {
    current = list.begin() + 1;
  }
  else current = list.begin();
}

////
void DeterministicIndices_2D::next_prune()
{
  ++current;
}

// remove current index and its pointer VERY TECHNIK for begin_r
void DeterministicIndices_2D::prune_current()
{
  current = list.erase(current);
}

////////////////
////////////////

unsigned int DeterministicIndices_2D::get_constraint_l()
{
  return *--current;
}

unsigned int DeterministicIndices_2D::get_constraint_r()
{
  return *--current;
}


////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// //// DeterministicIndices_2D2  ////
///// ///// ///// ///// ///// ///// /////

DeterministicIndices_2D2::DeterministicIndices_2D2() : Indices_2D2() {}

// full reset for pruning step
void DeterministicIndices_2D2::reset_prune()
{
  // reset iterators
  if (list.size() > 1)
  {
    current = list.begin() + 1;
  }
  else current = list.begin();
}


////
void DeterministicIndices_2D2::next_prune()
{
  ++current;
}

// remove current index and its pointer VERY TECHNIK for begin_r
void DeterministicIndices_2D2::prune_current()
{
  current = list.erase(current);
}

////////////////
////////////////

unsigned int DeterministicIndices_2D2::get_constraint_l()
{
  return *--current;
}

unsigned int DeterministicIndices_2D2::get_constraint_r()
{
  return list.back();
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// ///// RandomIndices_2D ///// /////
///// ///// ///// ///// ///// ///// /////

RandomIndices_2D::RandomIndices_2D() : Indices_2D2() {}

// full reset for pruning step,
void RandomIndices_2D::reset_prune()
{
  if (list.size() > 1){current = list.begin() + 1;}
  else current = list.begin();
}

// full next for pruning step
void RandomIndices_2D::next_prune() { ++current; }

// remove current index and its pointer
void RandomIndices_2D::prune_current() { current = list.erase(current); }



////////////////
////////////////

unsigned int RandomIndices_2D::get_constraint_l()
{
  unsigned int nbC_l = current - list.begin();
  return list[std::floor(dist(rng) * nbC_l)];
}


unsigned int RandomIndices_2D::get_constraint_r()
{
  unsigned int nbC_l = current - list.begin();
  return list[std::floor(dist(rng) * nbC_l)];
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// ///// RandomIndices_2D2 ///// /////
///// ///// ///// ///// ///// ///// /////

RandomIndices_2D2::RandomIndices_2D2() : Indices_2D2() {}

// full reset for pruning step,
void RandomIndices_2D2::reset_prune()
{
  if (list.size() > 1){current = list.begin() + 1;}
  else current = list.begin();
}

// full next for pruning step
void RandomIndices_2D2::next_prune() { ++current; }

// remove current index and its pointer
void RandomIndices_2D2::prune_current() { current = list.erase(current); }


////////////////
////////////////

unsigned int RandomIndices_2D2::get_constraint_l()
{
  unsigned int nbC_l = current - list.begin();
  return list[std::floor(dist(rng) * nbC_l)];
}

////////////////
////////////////

unsigned int RandomIndices_2D2::get_constraint_r()
{
  unsigned int nbC_r = list.end() - current - 1;
  return list[list.size() - std::ceil(dist(rng) * nbC_r)];
}




