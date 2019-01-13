  //Data struct we will use to perform the clustering

#ifndef NEWTYPES_H
#define NEWTYPES_H
using namespace std;

struct Laio
  {
    int snapshot;
    double rho;
    int cluster;
    int nnhd;
    double delta;
  };
  
    bool mycmp_rho(const Laio& a, const Laio& b)
  {
      return a.rho > b.rho;
  };

  
  //Data struct we will use to store the time, cv values, population and width of each cluster
  struct values {
  int iteration;
  int rank;
  double time;
  vector<double> cvs;
  int  population;
  vector<double> sigma;
  values(int iteration, int rank, double time,vector<double> cvs, int population, vector<double> sigma) 
    {
      this -> iteration = iteration;
      this -> rank = rank;
      this -> time = time;
      this -> cvs = cvs;
      this -> population = population;
      this -> sigma = sigma;
    }
  };
#endif	// NEWTYPES_H

