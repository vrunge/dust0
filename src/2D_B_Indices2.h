#ifndef INDICES2D2_H
#define INDICES2D2_H

#include <Rcpp.h>
#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;

class Indices_2D2
{
public:
  Indices_2D2();
  virtual ~Indices_2D2();

  // --- // Methods // --- //
  void reset();
  void next();
  bool check();

  void set_init_size(const unsigned int& size);
  void add(const unsigned int& value);

  std::vector<unsigned int> get_list();
  void remove_last();

  /////// --- // virtual Methods // --- ///////
  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;

  virtual unsigned int get_constraint_l() = 0;
  virtual unsigned int get_constraint_r() = 0;

  std::vector<unsigned int>::iterator current;

protected:
  std::vector<unsigned int> list;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

class DeterministicIndices_2D : public Indices_2D2
{

public:
  DeterministicIndices_2D();

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;

  unsigned int get_constraint_l() override;
  unsigned int get_constraint_r() override;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

class DeterministicIndices_2D2 : public Indices_2D2
{

public:
  DeterministicIndices_2D2();

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;

  unsigned int get_constraint_l() override;
  unsigned int get_constraint_r() override;

private:
  std::vector<unsigned int>::iterator begin_l;
  std::vector<unsigned int>::iterator begin_r;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

class RandomIndices_2D : public Indices_2D2
{

public:
  RandomIndices_2D();

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;

  unsigned int get_constraint_l() override;
  unsigned int get_constraint_r() override;

private:
  //////////// RANDOM NUMBER GENERATOR ////////////
  std::minstd_rand0 rng;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

class RandomIndices_2D2 : public Indices_2D2
{

public:
  RandomIndices_2D2();

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;

  unsigned int get_constraint_l() override;
  unsigned int get_constraint_r() override;

private:
  //////////// RANDOM NUMBER GENERATOR ////////////
  std::minstd_rand0 rng;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)
};


#endif
