#include <Rcpp.h>
#include <cmath>

#include "2D_B_Indices2.h"

using namespace Rcpp;

// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices_2D2::~Indices_2D2() {}

void Indices_2D2::reset(){current = list.begin();}
void Indices_2D2::next(){++current;}
bool Indices_2D2::check(){return current != list.end();}

unsigned int Indices_2D2::get_current(){return *current;}
std::forward_list<unsigned int> Indices_2D2::get_list(){return list;}
void Indices_2D2::remove_first(){list.pop_front();}

// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

RandomIndices_2D2::RandomIndices_2D2(unsigned int size, double alpha)
{

  // length of the random vector
  double k = std::max(2., ceil(pow(size, .2)));
  int len = std::log(alpha) / std::log(1 - 1/k);

  randomU = std::vector<double>(len, 0.);

  // Random number engine
  std::minstd_rand0  engine(std::random_device{}());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for(int i = 0; i < len; ++i) {randomU[i] = dist(engine);}
  u = randomU.begin();
}

void RandomIndices_2D2::add(unsigned int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void RandomIndices_2D2::reset_prune()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin();

  nbC = nb - 1;
}

void RandomIndices_2D2::next_prune()
{
  before = current;
  ++current;
  ++pointersCurrent;
  new_constraint();
}

void RandomIndices_2D2::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<unsigned int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
  new_constraint();
}

// --- // If no constraint can be selected, exit loop // --- //
bool RandomIndices_2D2::check_prune()
{
  return nbC > 0;
}

void RandomIndices_2D2::prune_last()
{
  list.erase_after(before);
  pointers.erase(std::next(pointersCurrent).base());
  nb--;
}

// --- // Select new constraint // --- //
// Optimisation possible, car dans le cas random on sélectionne une nouvelle contrainte avant de vérifier qu'elle sera utilisée
void RandomIndices_2D2::new_constraint()
{
  nbC--;
}

unsigned int RandomIndices_2D2::get_constraint()
{
  constraint = pointers[floor(nbC * (*u))];

  ++u;
  if (u == randomU.end())
  {
    u = randomU.begin();
  }

  return *constraint;
}


// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void DeterministicIndices_2D2::add(unsigned int value)
{
  list.push_front(value);
}

void DeterministicIndices_2D2::reset_prune()
{
  before = list.before_begin();
  current = std::next(before);
  constraint = std::next(current);
}

void DeterministicIndices_2D2::next_prune()
{
  before = current;
  current = constraint;
  new_constraint();
}

void DeterministicIndices_2D2::prune_current()
{
  current = list.erase_after(before);
  new_constraint();
}

bool DeterministicIndices_2D2::check_prune()
{
  return constraint != list.end();
}

void DeterministicIndices_2D2::prune_last()
{
  list.erase_after(before);
}

void DeterministicIndices_2D2::new_constraint()
{
  ++constraint;
}

unsigned int DeterministicIndices_2D2::get_constraint()
{
  return *constraint;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Random2Indices_2D2::Random2Indices_2D2(unsigned int size, double alpha)
{

  // length of the random vector
  double k = std::max(2., ceil(pow(size, .2)));
  int len = std::log(alpha) / std::log(1 - 1/k);

  randomU = std::vector<double>(len, 0.);

  // Random number engine
  std::minstd_rand0  engine(std::random_device{}());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for(int i = 0; i < len; ++i) {randomU[i] = dist(engine);}
  u = randomU.begin();
}

void Random2Indices_2D2::add(unsigned int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void Random2Indices_2D2::reset_prune()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin();

  nbC = nb - 1;
}

void Random2Indices_2D2::next_prune()
{
  before = current;
  ++current;
  ++pointersCurrent;
  new_constraint();
}

void Random2Indices_2D2::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<unsigned int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
  new_constraint();
}

// --- // If no constraint can be selected, exit loop // --- //
bool Random2Indices_2D2::check_prune()
{
  return nbC > 0;
}

void Random2Indices_2D2::prune_last()
{
  list.erase_after(before);
  pointers.erase(std::next(pointersCurrent).base());
  nb--;
}

// --- // Select new constraint // --- //
// Optimisation possible, car dans le cas random on sélectionne une nouvelle contrainte avant de vérifier qu'elle sera utilisée
void Random2Indices_2D2::new_constraint()
{
  nbC--;
}

unsigned int Random2Indices_2D2::get_constraint()
{
  constraint = pointers[floor(nbC * (*u))];

  ++u;
  if (u == randomU.end())
  {
    u = randomU.begin();
  }

  return *constraint;
}


// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void Deterministic2Indices_2D2::add(unsigned int value)
{
  list.push_front(value);
}

void Deterministic2Indices_2D2::reset_prune()
{
  before = list.before_begin();
  current = std::next(before);
  constraint = std::next(current);
}

void Deterministic2Indices_2D2::next_prune()
{
  before = current;
  current = constraint;
  new_constraint();
}

void Deterministic2Indices_2D2::prune_current()
{
  current = list.erase_after(before);
  new_constraint();
}

bool Deterministic2Indices_2D2::check_prune()
{
  return constraint != list.end();
}

void Deterministic2Indices_2D2::prune_last()
{
  list.erase_after(before);
}

void Deterministic2Indices_2D2::new_constraint()
{
  ++constraint;
}

unsigned int Deterministic2Indices_2D2::get_constraint()
{
  return *constraint;
}
