// [[Rcpp::depends(Rcpp)]]

#ifndef INDICES_H
#define INDICES_H

#include <Rcpp.h>
#include <forward_list>

using namespace Rcpp;

class IndicesModule
{
public:
  virtual ~IndicesModule();
  // --- // Methods // --- //
  virtual void add(int value) = 0;  // Pure virtual function

  void reset();
  void next();
  bool check();

  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual bool check_prune() = 0;

  virtual void prune_current() = 0;

  int get_current();
  virtual void new_constraint() = 0;
  virtual int get_constraint() = 0;
  std::forward_list<int> get_list();

  // --- // Fields // --- //
  std::forward_list<int> list;
  std::forward_list<int>::iterator current;
  std::forward_list<int>::iterator before;
};

class I_Random : public IndicesModule
{
public:
  I_Random(int size, double alpha = 1e-9);
  // ~I_Random() override;

  void add(int value) override;

  void reset_prune() override;
  void next_prune() override;
  bool check_prune() override;

  void prune_current() override;

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


class I_Deterministic : public IndicesModule
{
public:
  void add(int value) override;

  void reset_prune() override;
  void next_prune() override;
  bool check_prune() override;

  void prune_current() override;

  void new_constraint() override;
  int get_constraint() override;

private:
  std::forward_list<int>::iterator constraint;
};

#endif
