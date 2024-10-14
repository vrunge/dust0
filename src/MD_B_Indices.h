#ifndef INDICESMD_H
#define INDICESMD_H

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <xoshiro.h>

#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;


class Indices_MD
{

public:
  Indices_MD();
  Indices_MD(const unsigned int& nb_max_);
  virtual ~Indices_MD();

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
  virtual std::vector<unsigned int> get_constraints() = 0;

  std::vector<unsigned int> list;
  std::vector<unsigned int>::iterator current;

  void remove_last();

protected:
  unsigned int nb_max;

};


class VariableIndices_MD : public Indices_MD
{

public:
  VariableIndices_MD();
  VariableIndices_MD(const unsigned int& nb_max_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints() override;

private:
  std::vector<unsigned int>::iterator begin;

};


class RandomIndices_MD : public Indices_MD
{

public:
  RandomIndices_MD();
  RandomIndices_MD(const unsigned int& nb_max_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints() override;

private:
  dqrng::xoshiro256plus rng;
  std::uniform_real_distribution<double> dist;
  // std::uniform_int_distribution<int> dist2;
  std::vector<unsigned int>::iterator begin;

};

#endif
