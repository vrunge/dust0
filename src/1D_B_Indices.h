#ifndef INDICES1C_H
#define INDICES1C_H

#include <Rcpp.h>

#include <forward_list>

using namespace Rcpp;

class Indices
{
public:
  virtual ~Indices();
  
  // --- // Methods // --- //
  virtual void add(int value) = 0;
  
  void reset();
  void next();
  bool check();
  
  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;
  virtual bool check_prune() = 0;
  
  virtual void prune_last() = 0;
  
  int get_current();
  virtual void new_constraint() = 0;
  virtual int get_constraint() = 0;
  std::forward_list<int> get_list();
  
  // --- // Fields // --- //
  std::forward_list<int> list;
  std::forward_list<int>::iterator current;
  std::forward_list<int>::iterator before;
};

class RandomIndices : public Indices
{
public:
  RandomIndices(int size, double alpha = 1e-9);
  // ~I_Random() override;
  
  void add(int value) override;
  
  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;
  
  void prune_last() override;
  
  void new_constraint() override;
  int get_constraint() override;
  
private:
  int nb = 0;
  int nbC = 0;
  
  int* constraint;
  std::vector<int*> pointers;
  std::vector<int*>::reverse_iterator pointersCurrent;
  
  NumericVector randomU;
  NumericVector::iterator u;
};


class DeterministicIndices : public Indices
{
public:
  void add(int value) override;
  
  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;
  
  void prune_last() override;
  
  void new_constraint() override;
  int get_constraint() override;
  
private:
  std::forward_list<int>::iterator constraint;
};

#endif
