#ifndef RUN_STATS_HPP__
#define RUN_STATS_HPP__



class RunStats {
 protected:
  std::vector<double> vals;
  unsigned int n;

  double smean; // sample mean
  double svar; // sample variance
  double ssd; // sample standard deviation

  void update(const double & add);
  void update(const std::vector<double> & add);

 public:
  RunStats() : RunStats(std::vector<double>()) { };
  RunStats(const std::vector<double> & init);

  void operator () (const std::vector<double> & add){
    std::for_each(add.begin(),add.end(),
		  [this](const double & x){
		    update(x);
		  });
  };
				       
  void operator () (const double & add){
    update(add);
  };

  double mean() const {
    return mean;
  };

  double var() const {
    return var;
  };

  double sd() const {
    return std::sqrt(var);
  };

  
}




#endif
