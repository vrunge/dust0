// [[Rcpp::depends(RcppArmadillo)]]

#ifndef MC_H
#define MC_H

#include <RcppArmadillo.h>
#include <forward_list>
#include <vector>

using namespace Rcpp;

class MCHandler {
public:
  MCHandler();
  MCHandler(int constraints);
  ~MCHandler();
  void add(int value);
  
  void set_mu();
  
  void reset();
  void next();
  bool check();
  
  void reset_prune();
  void next_prune();
  bool check_prune();
  
  void reset_constraint();
  void next_constraint();
  bool check_constraint();
  
  bool allow();
  
  void prune();
  void simple_prune();
  
  int read();
  int read_constraint();
  double read_mu();
  double read_scale();
  int read_constraints();
  
  std::forward_list<int> get_list();
  
private:
  int nb = 0;
  int nbEnd;
  double scale;
  arma::colvec allMu;
  arma::colvec::iterator mu;
  
  std::forward_list<int> list;
  std::forward_list<int>::iterator current;
  std::forward_list<int>::iterator before;
  std::forward_list<int>::iterator constraintFrontier;
  std::forward_list<int>::iterator constraintIt;
};

#endif
