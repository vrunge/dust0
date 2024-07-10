// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "RandomList.h"
#include <forward_list>
#include <vector>

using namespace Rcpp;

// --------- // RandomList // --------- //
//
// This class allows drawing elements of a forwardlist at random, uniformly, among
// the elements following the current position of the forwardlist iterator created
// by RandomList::reset_prune() and incremented by RandomList::next_prune()
// This class was tailored to be used in the DUST algorithm
// 
// Parameters:


// Expliciter le constructeur, avec initialisation des 2 premiers indices
RandomList::RandomList(int inputSize, double alpha) {
  double k = std::max(2., ceil(pow(inputSize, .2)));
  randomU = Rcpp::runif(log(alpha) / log(1 - 1/k));
  u = randomU.begin();
  // reset_uniform(floor(pow(inputSize, 1.3)));
}

// Vérifier les fuites de mémoire
RandomList::~RandomList() {}

void RandomList::add(int value) {
  list.push_front(value);
  pointers.push_back(&list.front());
  length++;
}

void RandomList::next() {
  ++current;
}

void RandomList::next_prune() {
  before = current;
  ++current;
  ++pointersCurrent;
  lengthConstraints--;
}

void RandomList::reset() {
  current = list.begin();
}

void RandomList::reset_prune() {
  lengthConstraints = length - 1;
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin();
}

void RandomList::prune() {
  current = list.erase_after(before);
  pointersCurrent = std::vector<int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  length--;
  lengthConstraints--;
}

int RandomList::read() {
  return *current;
}

int RandomList::debug_read() {
  if (!check()) return -1;
  return *current;
}

int* RandomList::draw() {
  int* output = pointers[floor(lengthConstraints * (*u))];
  ++u;
  if (u == randomU.end())
  {
    // reset_uniform(floor(pow(inputSize, 1.3)));
    u = randomU.begin();
  }
  return output;
}

bool RandomList::check() {
  if (current == list.end()) {return false;}
  return true;
}

bool RandomList::check_prune() {
  // if (length <= 1 || lengthConstraints <= 0) {return false;}
  if (lengthConstraints <= 0) {return false;}
  return true;
}

std::forward_list<int> RandomList::get_list() {
  return list;
}


RCPP_MODULE(RandomListModule) {
  class_<RandomList>("ForwardListHandler")
  .constructor()
  .method("add", &RandomList::add)
  .method("next", &RandomList::next)
  .method("next_prune", &RandomList::next_prune)
  .method("reset", &RandomList::reset)
  .method("reset_prune", &RandomList::reset_prune)
  .method("check", &RandomList::check)
  .method("check_prune", &RandomList::check_prune)
  .method("prune", &RandomList::prune)
  .method("read", &RandomList::read)
  .method("debug_read", &RandomList::debug_read)
  .method("draw", &RandomList::draw)
  .method("get_list", &RandomList::get_list)
  ;
}
