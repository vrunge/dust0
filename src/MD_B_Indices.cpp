#include "MD_B_Indices.h"
#include <cmath>
using namespace Rcpp;

// --------- // Indices_MD // --------- //
//
// This class allows drawing elements of a forwardlist at random, uniformly, among
// the elements following the current position of the forwardlist iterator created
// by Indices_MD::reset_prune() and incremented by Indices_MD::next_prune()
// This class was tailored to be used in the DUST algorithm
//
// Parameters:

// Expliciter le constructeur, avec initialisation des 2 premiers indices
// Initiate handler and nb_max_ vector
Indices_MD::Indices_MD() : nb_l(1), nb_r(0), nb_max(1) {}
Indices_MD::Indices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_) :
  nb_l(nb_l_), nb_r(nb_r_), nb_max(nb_l_ + nb_r_)  {}

Indices_MD::~Indices_MD() {}

void Indices_MD::reset() { current = list.begin(); }
void Indices_MD::next() { ++current; }
bool Indices_MD::check() { return current != list.end(); }

void Indices_MD::set_init_size(const unsigned int& size) { list.reserve(size); }
void Indices_MD::add(const unsigned int& value){list.push_back(value);}

unsigned int Indices_MD::get_first() { return list.back(); }
std::vector<unsigned int> Indices_MD::get_list() { return list; }
void Indices_MD::remove_last() { list.pop_back(); }


////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// ////  DeterministicIndices_MD  ////
///// ///// ///// ///// ///// ///// /////

DeterministicIndices_MD::DeterministicIndices_MD() : Indices_MD() {}

DeterministicIndices_MD::DeterministicIndices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_) :
  Indices_MD(nb_l_, nb_r_)  { }

// full reset for pruning step
void DeterministicIndices_MD::reset_prune()
{
  // reset iterators
  if (list.size() > 1)
  {
    begin_l = list.begin();
    current = begin_l + 1;
    if (list.size() >= nb_r + 2) begin_r = list.end() - nb_r;
    else begin_r = current + 1;
  }
  else current = list.end();
}

////
void DeterministicIndices_MD::next_prune()
{
  ++current;
  /// move begin_l if too many indices bwn begin_l and current
  if (current - begin_l > nb_l) ++begin_l;
}

// remove current index and its pointer VERY TECHNIK for begin_r
void DeterministicIndices_MD::prune_current()
{
  unsigned int gap_r = begin_r - current;
  current = list.erase(current);
  begin_r = current + gap_r - 1;
}

////////////////
////////////////

std::vector<unsigned int> DeterministicIndices_MD::get_constraints_l()
{
  return std::vector(begin_l, current);
}

std::vector<unsigned int> DeterministicIndices_MD::get_constraints_r()
{
  if (begin_r == current) ++begin_r;
  return std::vector(begin_r, list.end());
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///// ///// ///// ///// ///// ///// /////
///// ///// RandomIndices_MD ///// /////
///// ///// ///// ///// ///// ///// /////

RandomIndices_MD::RandomIndices_MD() : Indices_MD() {}
RandomIndices_MD::RandomIndices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_)
  : Indices_MD(nb_l_, nb_r_)
  , rng(std::random_device{}())
  , dist(std::uniform_real_distribution(0., 1.))
{}



// full reset for pruning step,
void RandomIndices_MD::reset_prune()
{
  if (list.size() > 1){current = list.begin() + 1;} /// FAIRE SLT  current = list.begin() + 1;
  else current = list.begin();
}

// full next for pruning step
void RandomIndices_MD::next_prune() { ++current; }

// remove current index and its pointer
void RandomIndices_MD::prune_current() { current = list.erase(current); }



////////////////
////////////////

std::vector<unsigned int> RandomIndices_MD::get_constraints_l()
{
  std::vector<unsigned int> constraints;
  constraints.reserve(nb_l);
  unsigned int nbC_l = current - list.begin();
  for (unsigned int i = 0; i < nb_l; i++)
    constraints.push_back(list[std::floor(dist(rng) * nbC_l)]);

  return constraints;
}

////////////////
////////////////

std::vector<unsigned int> RandomIndices_MD::get_constraints_r()
{
  std::vector<unsigned int> constraints;
  constraints.reserve(nb_r);
  unsigned int nbC_r = list.end() - current - 1;
  if (nbC_r == 0) return constraints;
  for (unsigned int i = 0; i < nb_r; i++)
    constraints.push_back(list[list.size() - std::ceil(dist(rng) * nbC_r)]);
  return constraints;
}

//[[Rcpp::export]]
void test_indices()
{
  // // Initialize deterministic indices
  // DeterministicIndices_MD detIndices(2, 3);
  // detIndices.set_init_size(10);
  //
  // // Add values to the list
  // for (unsigned int i = 0; i < 1; ++i)
  // {
  //   detIndices.add(i + 1);
  // }
  //
  // // Test DeterministicIndices_MD
  // std::cout << "Testing DeterministicIndices_MD\n";
  // detIndices.reset_prune();
  // while (detIndices.check())
  // {
  //   std::cout << "Current: " << *detIndices.current << "\n";
  //
  //   std::vector<unsigned int> constraints_l = detIndices.get_constraints_l();
  //   std::vector<unsigned int> constraints_r = detIndices.get_constraints_r();
  //
  //   // std::cout << "Constraints L: ";
  //   // for (const auto& val : constraints_l) std::cout << val << " ";
  //   // std::cout << "\n";
  //   //
  //   // std::cout << "Constraints R: ";
  //   // for (const auto& val : constraints_r) std::cout << val << " ";
  //   // std::cout << "\n";
  //   //
  //   // detIndices.next_prune();
  // }

  // Initialize random indices
  RandomIndices_MD randIndices(2, 3);
  randIndices.set_init_size(10);

  // Add values to the list
  for (unsigned int i = 0; i < 1; ++i)
  {
    randIndices.add(i + 1);
  }

  // Test RandomIndices_MD
  std::cout << "\nTesting RandomIndices_MD\n";
  randIndices.reset_prune();
  while (randIndices.check())
  {
    std::cout << "Current: " << *randIndices.current << "\n";

    std::vector<unsigned int> constraints_l = randIndices.get_constraints_l();
    std::vector<unsigned int> constraints_r = randIndices.get_constraints_r();

    std::cout << "Constraints L: ";
    for (const auto& val : constraints_l) std::cout << val << " ";
    std::cout << "\n";

    std::cout << "Constraints R: ";
    for (const auto& val : constraints_r) std::cout << val << " ";
    std::cout << "\n";

    randIndices.next_prune();
  }

  std::vector<unsigned> test = { 0 };
  std::vector<unsigned> test2(test.begin() + 1, test.end());
}

