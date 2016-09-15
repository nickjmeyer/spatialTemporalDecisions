#include <gflags/gflags.h>
#include <glog/logging.h>
#include "system.hpp"
#include "model.hpp"
#include "model2GravityEDist.hpp"
#include "model2EdgeToEdge.hpp"

DEFINE_string(srcDir,"","Path to source directory");
DEFINE_bool(edgeToEdge,false,"Edge to edge transmission");
DEFINE_string(outDir,"","Path to save parameters");


template <typename T>
void copyParams(const boost::filesystem::path path) {
	System<T,T> s;
	s.modelGen_r.read();
	s.modelGen_r.save_to(path);
}


int main(int argc, char ** argv) {
  ::google::InitGoogleLogging(argv[0]);
  ::google::ParseCommandLineFlags(&argc,&argv,true);
  njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

	const boost::filesystem::path path(FLAGS_outDir);

	if(FLAGS_edgeToEdge) {
		copyParams<Model2EdgeToEdge>(path);
	} else {
    copyParams<Model2GravityEDist>(path);
	}
}
