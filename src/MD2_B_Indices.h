#ifndef INDICESMD2_H
#define INDICESMD2_H

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <xoshiro.h>

#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;


class Indices_MD2
{

public:
  Indices_MD2();
  Indices_MD2(const unsigned int& nb_l_, const unsigned int& nb_r_);
  virtual ~Indices_MD2();

  // --- // Methods // --- //
  void reset();
  void next();
  bool check();

  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;
  virtual bool check_prune() = 0;

  void add(const unsigned int& value);

  void set_size(const unsigned int& size);
  std::vector<unsigned int> get_list();
  virtual std::vector<unsigned int> get_constraints_l() = 0;
  virtual std::vector<unsigned int> get_constraints_r() = 0;

  /// WHY NOT IN PROTECTED?
  std::vector<unsigned int> list;
  std::vector<unsigned int>::iterator current;

  void remove_last();

protected:
  unsigned int nb_l;
  unsigned int nb_r;
  unsigned int nb_max;

};


////////////////////////////////////////////////
////////////////////////////////////////////////

class VariableIndices_MD2 : public Indices_MD2
{

public:
  VariableIndices_MD2();
  VariableIndices_MD2(const unsigned int& nb_l_, const unsigned int& nb_r_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints_l() override;
  std::vector<unsigned int> get_constraints_r() override;

private:
  std::vector<unsigned int>::iterator begin;

};

////////////////////////////////////////////////
////////////////////////////////////////////////

class RandomIndices_MD2 : public Indices_MD2
{

public:
  RandomIndices_MD2();
  RandomIndices_MD2(const unsigned int& nb_l_, const unsigned int& nb_r_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints_l() override;
  std::vector<unsigned int> get_constraints_r() override;

private:
  dqrng::xoshiro256plus rng;
  std::uniform_real_distribution<double> dist;
  // std::uniform_int_distribution<int> dist2;
  std::vector<unsigned int>::iterator begin;

};

#endif
