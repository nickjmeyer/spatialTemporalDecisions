#ifndef SETTINGS_HPP__
#define SETTINGS_HPP__

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include "utilities.hpp"
#include "dataDepth.hpp"

class Settings {
 public:
  Settings();
  ~Settings();

  void timeElapsed();
  
  void set(int numInitVals, char ** initVals);

  void clean();
  int cleaned;

  std::string fileName;
  std::string srcDir;
  std::string datDir;
  std::string date;

  const static int numVals;
  
  int tick,tock,seconds;
  void timeStamp();

  std::string datExt(int inDir, std::string beg, std::string end) const;
  std::string datExt(std::string beg, std::string end) const;
  std::string datExt(int inDir, std::string end) const;
  std::string datExt(std::string end) const;
  std::string datExt(int inDir) const;
  std::string datExt() const;

  std::string srcExt(const std::string file);
};

namespace njm{
  extern Settings sett;
};

#endif
