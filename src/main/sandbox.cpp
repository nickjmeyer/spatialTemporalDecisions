#include <gflags/gflags.h>
#include <glog/logging.h>
#include "runM1Mles.hpp"

using namespace google;
using namespace gflags;


int main(int argc, char ** argv){
    InitGoogleLogging(argv[0]);
    ParseCommandLineFlags(&argc,&argv,true);

    std::vector<int> fips;
    njm::fromFile(fips,"data/wns/fips.txt");
    const int numNodes = fips.size();
    CHECK_EQ(1128,numNodes);

    std::vector< std::vector<int> > obsData;
    {
        std::vector<int> obsData_raw;
        njm::fromFile(obsData_raw,"data/wns/obsData.txt");

        std::vector<int> add;
        std::vector<int>::const_iterator it,end;
        it = obsData_raw.begin();
        end = obsData_raw.end();
        int count = 0;
        while(it != end) {
            add.push_back(*it);
            ++it;
            ++count;
            if (count == numNodes) {
                obsData.push_back(add);
                add.clear();
                count = 0;
            }
        }
        CHECK_EQ(0,count);
    }

    for (int i = 0; i < (obsData.size() - 1); ++i) {
        for (int j = 0; j < numNodes; ++j) {
            const int before = obsData.at(i).at(j)/2;
            const int after = obsData.at(i+1).at(j)/2;
            CHECK_LE(before,after)
                << "node " << j << " failed for time point "
                << i << " to " << i+1;
        }
    }

    return 0;
}
