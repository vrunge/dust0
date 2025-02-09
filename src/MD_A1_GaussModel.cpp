#include "MD_A1_GaussModel.h"

using namespace Rcpp;

Gauss_MD::Gauss_MD(int dual_max_type, int constraints_type, Nullable<unsigned> nbLoops)
  : DUST_MD(dual_max_type, constraints_type, nbLoops) {}

double Gauss_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double result = 0;
  for (unsigned int row = 0; row < dim; row++)
    result += pow(cumsum(row, t) - cumsum(row, s), 2);

  return - .5 * result / (t - s);
}

double Gauss_MD::statistic(const double& value) const
{
  return value;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// GAUSS function =>  -0.5 sum (a_i + X b_i)^2/(c + X d) -e X - f

double Gauss_MD::dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  double result = 0;
  for (unsigned int i = 0; i <  dim; i++){result += pow(a(i) + point*b(i), 2);}
  result = - 0.5 * pow(c + point * d, -1) * result;
  return(result - e*point - f);
}

////////////////////////////////////////////////////////////////////////////////

//// SEARCH FOR THE BEST value in one direction on the DUAL FUNCTION
//
// GAUSS function =>  -0.5 sum (a_i + X b_i)^2/(c + X D) -e X - f
// write D to be different from dimension d
//
double Gauss_MD::dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  double Max = -std::numeric_limits<double>::infinity();

  /// Find mu_min - mu_max
  std::array<double, 2>  mu_interval = muInterval(a,b,c,d);

  // GAUSS function =>  -0.5 sum (a_i + X b_i)^2/(c + X d) -e X - f
  // discriminant of the 2nd order polynomial (after differentiating):
  // roots of A mu^2 + B mu + C ?
  // = delta * DELTA2
  //
  double A2 = arma::dot(a,a);
  double B2 = arma::dot(b,b);
  double AB = arma::dot(a,b);

  double DELTA2 = d*d*A2 - 2*c*d*AB + c*c*B2;
  double delta = B2 + 2*d*e;

  if(delta < 0)  /// No root case
  {
    if(d > 0)    ///  A mu^2 + B mu + C > 0, max in +INFTY
    {
      //Rcout << "D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 D > 0 " << std::endl;
      argmax = mu_interval[1];
      //Rcout << argmax << std::endl;
      //Rcout << Max << std::endl;
      //Rcout << D << std::endl;
      //Rcout << delta << std::endl;
      //Rcout << DELTA2 << std::endl;
      //throw std::runtime_error("Something went wrong! Stopping execution.");
      return(std::numeric_limits<double>::infinity());
    }
    else if(d < 0)  ///  A mu^2 + B mu + C < 0, max in 0
    {
      //Rcout << "D < 0 D < 0 D < 0 D < 0 D < 0 D < 0 " << std::endl;
      //Rcout << "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD" << std::endl;
      //Rcout << c << " / " << D << " / " << e << " / " << f << std::endl;
      argmax = mu_interval[0];
      Max = dual1D_Eval(argmax,a,b,c,d,e,f);
      //Rcout << Max << std::endl;
      //Rcout << argmax << std::endl;
      //Rcout << Max << std::endl;
      return(Max);
    }
  }
  else
  {
    //Rcout << "root root root root root root root root root " << std::endl;
    double root = (-c + std::sqrt(DELTA2/delta))/d;
    argmax = (0 > root) ? 0 : root;
    Max = dual1D_Eval(argmax,a,b,c,d,e,f);
    //Rcout << Max << std::endl;
    return(Max);
  }
  return(Max);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_MD::muMax(const double& a, const double& b) const
{
  return 1.;
}


std::array<double, 2> Gauss_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const
{
  std::array<double, 2> interval = {0, std::numeric_limits<double>::infinity()};
  if (c > 0 && d < 0)
  {
    interval[1] = -c / d;
  }
  else if (c < 0 && d > 0)
  {
    interval[0] = -c / d;
  }
  return(interval);
}



void Gauss_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  // DIRECTION_SUM MUST BE POSITIVE HERE
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Gauss_MD::Dstar(const double& x) const
{
  return .5 * pow(x, 2);
}

double Gauss_MD::DstarPrime(const double& x) const
{
  return x;
}

double Gauss_MD::DstarSecond(const double& x) const
{
  return 1.;
}

std::string Gauss_MD::get_model() const { return "gauss"; }





