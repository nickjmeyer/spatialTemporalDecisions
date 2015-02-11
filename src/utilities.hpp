#ifndef UTILITIES_HPP__
#define UTILITIES_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include <cmath>

namespace njm{
  
  template<class T>
  std::string toString(const T n, 
		       const std::string end = "",
		       const int w = 16, const int p = 8){
    std::stringstream ss;
    ss << std::setw(w) << std::setprecision(p) << std::fixed
       << n << end;
    return ss.str();
  };


  template<class T>
  std::string toString(const std::vector<T> & n, 
		       const std::string sep = "  ",
		       const std::string end = "",
		       const int w = 16, const int p = 8){
    std::stringstream ss;
    typename
    std::vector<T>::const_iterator it;
    for(it = n.begin(); it!=n.end(); it++)
      ss << std::setw(w) << std::setprecision(p) << std::fixed
	 << njm::toString(*it) << sep;
    ss << end;
    return ss.str();
  };


  template<class T>
  int toFile(const T & n,
	     const std::string file,
	     const std::ios_base::openmode mode = std::ios_base::app,
	     const std::string end = "\n",
	     const int w = 16, const int p = 8){
    std::ofstream ofs;
    ofs.open(file.c_str(), mode);
    if(ofs.is_open()){
      ofs << toString(n,end,w,p);
      ofs.close();
      return 0;
    }
    else
      return 1;
  };
  template<class T>
  int toFile(const std::vector<T> & n,
	     const std::string file,
	     const std::ios_base::openmode mode = std::ios_base::app,
	     const std::string sep = "  ",
	     const std::string end = "\n",
	     const int w = 16, const int p = 8){
    std::ofstream ofs;
    ofs.open(file.c_str(), mode);
    if(ofs.is_open()){
      ofs << toString(n,sep,end,w,p);
      ofs.close();
      return 0;
    }
    else
      return 1;
  };

  
  template<class T>
  int fromFile(T & n,
	       const std::string file){
    std::ifstream ifs;
    ifs.open(file.c_str(),std::ifstream::in);
    if(ifs.is_open()){
      ifs >> n;
      ifs.close();
      return 0;
    }
    else{
      std::cout << "Failed to read " << file << std::endl;
      return 1;
    }
  };


  template<class T>
  int fromFile(std::vector<T> & n,
	       const std::string file){
    std::ifstream ifs;
    ifs.open(file.c_str(),std::ifstream::in);
    if(ifs.is_open()){
      T t;
      n.clear();
      ifs >> t;
      while(!ifs.eof()){
	n.push_back(t);
	ifs >> t;
      }
      ifs.close();
      return 0;
    }
    else{
      std::cout << "Failed to read " << file << std::endl;
      return 1;
    }
  };


  template<class T>
  void message(const T & s){
    if(omp_get_thread_num() == 0)
      std::cout << '#' + std::string(63,'-') << std::endl
		<< "Message: " << std::endl
		<< njm::toString(s)
		<< std::endl
		<< std::string(63,'-') + '#' << std::endl;
  };

  double expit(const double n);


  double l2norm(const std::vector<double> & v0,
		const std::vector<double> & v1);

  double l2norm(std::vector<double> & v);

};

#endif
