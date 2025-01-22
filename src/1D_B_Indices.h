#ifndef INDICES1C_H
#define INDICES1C_H

#include <Rcpp.h>

#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;

class Indices_1D
{
public:
  virtual ~Indices_1D();

  // --- // Methods // --- //
  virtual void add(unsigned int value) = 0;

  void reset(); //iterator current = start
  void next(); //current++
  bool check(); //no the end

  unsigned int get_first();
  void remove_first(); //remove the first element
  std::forward_list<unsigned int> get_list();

  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;
  virtual bool check_prune() = 0;

  virtual void prune_last() = 0;

  unsigned int get_current();
  virtual void new_constraint() = 0;
  virtual unsigned int get_constraint() = 0;

protected:
  // --- // Fields // --- //
  std::forward_list<unsigned int> list;
  std::forward_list<unsigned int>::iterator current;
  std::forward_list<unsigned int>::iterator before;
};



////////////////
//////////////// derived class 1:  RandomIndices_1D
////////////////

class RandomIndices_1D : public Indices_1D
{
public:
  RandomIndices_1D();

  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  unsigned int nb = 0; // number of elements in the list = length
  unsigned int nbC = 0; // number of available constraint (smaller indices) for a given s index

  unsigned int* constraint; // the constraint r (for s, r < s)
  std::vector<unsigned int*> pointers; // all pointers for the list of indices
  std::vector<unsigned int*>::reverse_iterator pointersCurrent; // to move on pointers

  //////////// RANDOM NUMBER GENERATOR ////////////
  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)
};

////////////////
//////////////// derived class 1:  DeterministicIndices_1D
////////////////

class DeterministicIndices_1D : public Indices_1D
{
public:
  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  // iterator for list
  // a kind of "after" or "next" pointer
  std::forward_list<unsigned int>::iterator constraint;
};

#endif
