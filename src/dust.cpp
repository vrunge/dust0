#include <forward_list>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dust(NumericVector data, double penalty = 0)
{
  int n = data.size();
  if(penalty == 0){penalty = 2*log(n);}
  ///
  /// DATA TRANSFORMATION
  ///
  double cumsums[n + 1];
  cumsums[0] = 0;
  for (int i = 0; i < n; i++)
  {
    cumsums[i + 1] = cumsums[i] + data[i];
  }

  //for (int i = 0; i < n + 1; i++)
  //{
  //  Rcout << cumsums[i] << " --- ";
  //}
  //Rcout << std::endl;
  //Rcout << penalty << std::endl;


  ///
  /// INITIALISATION
  ///
  double costQ[n + 1]; //costQ[i] optimal cost for data y(1) to y(i-1)
  int cp[n + 1]; //cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  costQ[0] = -penalty;
  cp[0] = 0;

  ///
  std::forward_list<int> indexSet = {};
  int nb[n + 1]; //nb element in indexSet
  nb[0] = 0;

  std::forward_list<int>::iterator itr;
  for(itr=indexSet.begin();itr!=indexSet.end();++itr) Rcout <<*itr<<' ';
  Rcout <<'\n';

  std::forward_list<int> list = {1, 2, 3, 4, 5,1,2,3, 30000};
  int *arr = new int[9];
  std::copy(std::begin(list), std::end(list), arr);

  // Print the copied elements of the array
  for (int i = 0; i < 9; i++) {
    std::cout << arr[i] << " ";
  }
  delete []arr;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  ///
  /// ALGORITHM OP
  ///
  for (int i = 1; i < n + 1; i++) // using cumsums[i]
  {


    ///
    /// pruning
    ///
    int *my_copy = new int[nb[i]];
    std::copy(std::begin(indexSet), std::end(indexSet), my_copy);


    delete []my_copy;
  }

  return(data);
}

