// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "costFunctions.h"
#include "multipleConstraints.h"
#include <forward_list>
#include <cmath>
#include <random>
#include <limits>
using namespace Rcpp;


// --------- // Sum // --------- //
//
// Sums the values of input vector
//
// Parameters:
//  - data (vector): a vector of numeric values

// [[Rcpp::export]]
arma::colvec colPow(arma::colvec& inputCol, int power)
{
  arma::colvec out(inputCol.n_elem);
  arma::colvec::iterator iOut = out.begin();
  auto iIn = inputCol.begin();
  do
  {
    *iOut = pow(*iIn, power);
    ++iOut;
    ++iIn;
  }
  while (iOut != out.end());
  return out;
}



// --------- // Simple test // --------- //
// --------- // Simple test // --------- //
// --------- // Simple test // --------- //
//
// Performs the simple lagrangian duality test to find the lowest "visible " point
// of a given cost function
//
// Parameters:
//  - t: the current step in the OP main loop
//  - i: the current step in the pruning sub loop
//  - j: a random index drawn from the available indices, less than i
//  - valuesCumsum (vector): the cumulative arma::sum of the data (index 0 == 0)
//  - costRecord (vector): the optimal model cost at each OP index


// [[Rcpp::export]]
double test
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord
)
{

  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(i+1):t]; y_ji = y[(j+1):i]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi + mu * (Qi-Qj) - ((t - i) * m_it - mu * (i - j) * sm_ji)^2 / (t - i - mu * (i - j))
  // Denoting mu_max = (t - i) / (i - j)
  //
  // Formula: mu* = max(0, mu_max * (1 - abs(m_it - m_ji) / sqrt((Qi - Qj) / (i-j) + m_ji^2)))

  int objectiveLength = t - i;
  double objectiveSum = valuesCumsum[t] - valuesCumsum[i];
  double objectiveMean = objectiveSum / objectiveLength; // m_it
  double squareObjective = pow(objectiveMean, 2);

  int gapLength = i - j;
  double gapSum = valuesCumsum[i] - valuesCumsum[j];
  double gapMean = gapSum / gapLength; // m_ji
  double crossMean = 2 * objectiveMean * gapMean;

  double costTerm = (costRecord[i] - costRecord[j]) / gapLength; // Qi - Qj / i - j

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costTerm <= squareObjective - crossMean)
    return costRecord[i] - objectiveLength * squareObjective; // D(0)

  double squareGap = pow(gapMean, 2);
  double gapTerm = costTerm + squareGap;

  return
  costRecord[i]
  + objectiveLength
  * (
      gapTerm
      + pow(gapMean, 2)
    - crossMean
    - 2 * fabs(objectiveMean - gapMean) * sqrt(gapTerm)
  );
}


// [[Rcpp::export]]
double test1D
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
)
{

  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(i+1):t]; y_ji = y[(j+1):i]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi + (t - i) ( mu (Qi-Qj)/(i-j) - (m_it - mu m_ji)/(1-mu))
  //
  // Formula: mu* = max(0, 1 - abs(m_it - m_ji)/sqrt((Qi-Qj)/(i-j) + m_ji^2))
  // Formula: mu* > 0, d* = D(mu*) = Qi - (t-i) m_it^2 + (t-i) (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Formula: (d* - Qt) / (t - i) = (Qi-Qt)/(t-i) - m_it^2 + (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Pruning happens if (d* - Qt) / (t - i) > 0

  int objectiveLength = t - i;
  double objectiveMean = (valuesCumsum[t] - valuesCumsum[i]) / objectiveLength; // m_it

  int gapLength = i - j;
  double gapMean = (valuesCumsum[i] - valuesCumsum[j]) / gapLength; // m_ji

  double costJI = (costRecord[i] - costRecord[j]) / gapLength + pow(gapMean, 2); // Qi - Qj / i - j
  double costIT = (optimalCost - costRecord[i]) / objectiveLength + pow(objectiveMean, 2);
  double meanGap = fabs(gapMean - objectiveMean);

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costJI <= pow(meanGap, 2))
    return - costIT;

  // Case 2: mu* > 0
  return - costIT + pow(sqrt(costJI) - meanGap, 2);
}


// [[Rcpp::export]]
double test1DRand
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
)
{

  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(i+1):t]; y_ji = y[(j+1):i]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi - (t - i) m_it^2 + (t - i) ( mu ((Qi-Qj)/(i-j) + m_ji^2) - (m_it - mu m_ji)/(1-mu))
  //
  // test value is given by formula:
  //    (D(mu) - Qt) / (t-i) = (Qi - Qt) / (t - i) - m_it^2 + mu ((Qi - Qj) / (i - j) + m_ji^2) + (m_it - mu m_ji)^2 / (1 - mu)

  int objectiveLength = t - i;
  double objectiveMean = (valuesCumsum[t] - valuesCumsum[i]) / objectiveLength; // m_it

  int gapLength = i - j;
  double gapMean = (valuesCumsum[i] - valuesCumsum[j]) / gapLength; // m_ji

  double mu = runif(1)[0];

  return (costRecord[i] - optimalCost) / objectiveLength - pow(objectiveMean, 2)
    + mu * ((costRecord[i] - costRecord[j]) / gapLength + pow(gapMean, 2))
    - pow(gapMean - mu * objectiveMean, 2) / (1 - mu);
}


// [[Rcpp::export]]
double testMV1C
(
    int t,
    int i,
    int j,
    const arma::dmat& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
)
{

  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(i+1):t]; y_ji = y[(j+1):i]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi + (t - i) ( mu (Qi-Qj)/(i-j) - (m_it - mu m_ji)/(1-mu))
  //
  // Formula: mu* = max(0, 1 - abs(m_it - m_ji)/sqrt((Qi-Qj)/(i-j) + m_ji^2))
  // Formula: mu* > 0, d* = D(mu*) = Qi - (t-i) m_it^2 + (t-i) (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Formula: (d* - Qt) / (t - i) = (Qi-Qt)/(t-i) - m_it^2 + (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Pruning happens if (d* - Qt) / (t - i) > 0


  int objectiveLength = t - i;
  arma::colvec objectiveMean = (valuesCumsum.col(t) - valuesCumsum.col(i)) / objectiveLength; // m_it

  int gapLength = i - j;
  arma::colvec gapMean = (valuesCumsum.col(i) - valuesCumsum.col(j)) / gapLength; // m_ji

  double costJI = (costRecord[i] - costRecord[j]) / gapLength
  + arma::sum(colPow(gapMean));
  double costIT = (optimalCost - costRecord[i]) / objectiveLength
  + arma::sum(colPow(objectiveMean));

  arma::colvec meanGap = gapMean - objectiveMean;
  double sqMeanGap = arma::sum(colPow(meanGap));

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costJI <= sqMeanGap)
    return - costIT;

  // Case 2: mu* > 0
  return - costIT + pow(sqrt(costJI) - sqrt(sqMeanGap), 2);
}


double testMV
(
    int t,
    int i,
    MCHandler& indices,
    const arma::dmat& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
)
{
  int d = valuesCumsum.n_rows;

  arma::subview_col sumi = valuesCumsum.col(i);
  double costi = costRecord[i];

  int objectiveLength = t - i;
  arma::colvec objectiveMean = (valuesCumsum.col(t) - sumi) / objectiveLength; // m_it
  double costIT = (optimalCost - costRecord[i]) / objectiveLength
    + arma::sum(colPow(objectiveMean, 2));

  double constraintMean;

  double costTerm = 0.;
  double gapTerm = 0.;
  double squareGap;

  IntegerVector lengths(indices.read_constraints());
  IntegerVector::iterator lengthj;

  int j;
  double mu;
  double sumdi;
  for(int dim = 0; dim < d; dim++)
  {
    squareGap = objectiveMean(dim);
    sumdi = sumi(dim);
    lengthj = lengths.begin();
    indices.reset_constraint();
    do
    {
      j = indices.read_constraint();
      mu = indices.read_mu();

      *lengthj = i - j;
      costTerm += mu * (costi - costRecord[j]) / *lengthj;

      ++lengthj;
      indices.next_constraint();
    } while (indices.check_constraint());

    lengthj = lengths.begin();
    indices.reset_constraint();
    do
    {
      j = indices.read_constraint();
      mu = indices.read_mu();

      constraintMean = (sumdi - valuesCumsum(dim, j)) / *lengthj;

      costTerm -= mu * (pow(constraintMean, 2));
      squareGap -= mu * constraintMean;
      indices.next_constraint();
    } while (indices.check_constraint());
    gapTerm += pow(squareGap, 2);
  }

  return - costIT + costTerm - gapTerm / (1 - indices.read_scale());
}

// --------- // modelCost // --------- //
// --------- // modelCost // --------- //
// --------- // modelCost // --------- //
//
// Computes the gaussian model cost in the OP algorithm at step t with cursor at i
//
// Parameters:
//  - t: the current step in the OP main loop
//  - i: the current step in the OP sub loop
//  - valuesCumsum (vector): the cumulative arma::sum of the data (index 0 == 0)
//  - costRecord (vector): the optimal model cost at each OP index

// [[Rcpp::export]]
double modelCost
(
    int t,
    int i,
    NumericVector& valuesCumsum
)
{
  return - pow(valuesCumsum[t] - valuesCumsum[i], 2) / (t - i);
}

// [[Rcpp::export]]
double modelCostMV
(
    int t,
    int i,
    arma::dmat& valuesCumsum
)
{
  arma::colvec midResult = valuesCumsum.col(t) - valuesCumsum.col(i);
  return - arma::sum(colPow(midResult)) / (t - i);
}


// --------- // MAD Estimator // --------- //
// --------- // MAD Estimator // --------- //
// --------- // MAD Estimator // --------- //
//
// Provided some data vector, computes the Median Absolute Deviation estimator
//
// Parameters:
//  - data (vector): a vector of numeric values

// [[Rcpp::export]]
double madEstimator(NumericVector& data)
{
  int n = data.size();
  int medianIndex = n / 2 + n % 2;

  NumericVector sorted_data(data.begin(), data.end());
  std::sort(sorted_data.begin(), sorted_data.end());

  double median = sorted_data[medianIndex];

  double mad = 0.;
  NumericVector::iterator y = data.begin();

  do
  {
    mad += fabs(*y - median);
    ++y;
  }
  while (y != data.end());

  return mad / n;
}
