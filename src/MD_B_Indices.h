#ifndef INDICESMD_H
#define INDICESMD_H

#include <RcppArmadillo.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;

class Indices_MD
{

public:
  Indices_MD();
  Indices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_);
  virtual ~Indices_MD();

  // --- // Methods // --- //
  void reset();
  void next();
  bool check();

  void add(const unsigned int& value);
  void set_size(const unsigned int& size);

  std::vector<unsigned int> get_list();
  void remove_last();

  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;
  virtual bool check_prune() = 0;

  virtual std::vector<unsigned int> get_constraints_l() = 0;
  virtual std::vector<unsigned int> get_constraints_r() = 0;

  std::vector<unsigned int>::iterator current;


protected:
  std::vector<unsigned int> list;
  unsigned int nb_l;
  unsigned int nb_r;
  unsigned int nb_max;

};


////////////////////////////////////////////////
////////////////////////////////////////////////

class DeterministicIndices_MD : public Indices_MD
{

public:
  DeterministicIndices_MD();
  DeterministicIndices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints_l() override;
  std::vector<unsigned int> get_constraints_r() override;

private:
  std::vector<unsigned int>::iterator begin_l;
  std::vector<unsigned int>::iterator begin_r;
};

////////////////////////////////////////////////
////////////////////////////////////////////////

class RandomIndices_MD : public Indices_MD
{

public:
  RandomIndices_MD();
  RandomIndices_MD(const unsigned int& nb_l_, const unsigned int& nb_r_);

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  std::vector<unsigned int> get_constraints_l() override;
  std::vector<unsigned int> get_constraints_r() override;

private:
  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 rng;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)
};

#endif
