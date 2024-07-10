// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "multipleConstraints.h"
#include "logging.h"
#include <forward_list>
#include <vector>

using namespace Rcpp;

// --------- // MCHandler // --------- //
//
// This class allows drawing elements of a forwardlist at random, uniformly, among
// the elements following the current position of the forwardlist iterator created
// by MCHandler::reset_prune() and incremented by MCHandler::next_prune()
// This class was tailored to be used in the DUST algorithm
// 
// Parameters:


// Expliciter le constructeur, avec initialisation des 2 premiers indices
// Initiate handler and constraints vector
MCHandler::MCHandler() : nbEnd(1) {}
MCHandler::MCHandler(int constraints) : nbEnd(constraints) {}

// Vérifier les fuites de mémoire
MCHandler::~MCHandler() {}

// add index
void MCHandler::add(int value) {
  list.push_front(value);
  nb++;
}

// set mu values and scale them with random scale factor
void MCHandler::set_mu() {
  allMu = arma::randu<arma::colvec>(nbEnd);
  scale = runif(1)[0];
  allMu = scale * allMu / arma::sum(allMu);
}

// simple reset for OP step
void MCHandler::reset() {
  current = list.begin();
}

// simple next for OP step
void MCHandler::next() {
  ++current;
}

// simple break check for OP step
bool MCHandler::check() {
  return current != list.end();
}

// full reset for pruning step, including reinitialization of constraints vector
void MCHandler::reset_prune() {
  // reset iterators
  current = list.begin();
  before = list.before_begin();
  constraintFrontier = list.begin();
  std::advance(constraintFrontier, nbEnd);
  set_mu();
}

// full next for pruning step
void MCHandler::next_prune() {
  before = current;
  set_mu();
  ++current;
  ++constraintFrontier;
}

// remove current index and its pointer
void MCHandler::prune() {
  current = list.erase_after(before);
  set_mu();
  ++constraintFrontier;
  nb--;
}

// break check for pruning step
bool MCHandler::check_prune() {
  return constraintFrontier != list.end();
}

// reset constraint iterator
void MCHandler::reset_constraint()
{
  mu = allMu.begin();
  constraintIt = std::next(current);
}

// increment constraint iterator
void MCHandler::next_constraint()
{
  ++mu;
  ++constraintIt;
}

// break check for constraint loop
bool MCHandler::check_constraint()
{
  return mu != allMu.end();
}

// check if pruning step is possible
bool MCHandler::allow()
{
  return nb > nbEnd;
}

// prune without incrementing frontier
void MCHandler::simple_prune()
{
  current = list.erase_after(before);
  nb--;
}

// obtain current tested index
int MCHandler::read() {
  return *current;
}

// obtain constraint index
int MCHandler::read_constraint() {
  return *constraintIt;
}

// obtain current mu
double MCHandler::read_mu()
{
  return *mu;
}

// obtain the sum of mu values
double MCHandler::read_scale()
{
  return scale;
}

// obtain number of constraints to be considered
int MCHandler::read_constraints()
{
  return nbEnd;
}

// obtain list for main function output
std::forward_list<int> MCHandler::get_list() {
  return list;
}


RCPP_MODULE(MCHandlerModule) {
  class_<MCHandler>("MCHandler")
  .constructor()
  .constructor<int>()
  .method("add", &MCHandler::add)
  .method("set_mu", &MCHandler::set_mu)
  .method("reset", &MCHandler::reset)
  .method("next", &MCHandler::next)
  .method("check", &MCHandler::check)
  .method("reset_prune", &MCHandler::reset_prune)
  .method("next_prune", &MCHandler::next_prune)
  .method("check_prune", &MCHandler::check_prune)
  .method("reset_constraint", &MCHandler::reset_constraint)
  .method("next_constraint", &MCHandler::next_constraint)
  .method("check_constraint", &MCHandler::check_constraint)
  .method("prune", &MCHandler::prune)
  .method("simple_prune", &MCHandler::simple_prune)
  .method("read", &MCHandler::read)
  .method("read_constraint", &MCHandler::read_constraint)
  .method("read_mu", &MCHandler::read_mu)
  .method("read_scale", &MCHandler::read_scale)
  .method("read_constraints", &MCHandler::read_constraints)
  .method("get_list", &MCHandler::get_list)
  ;
}
