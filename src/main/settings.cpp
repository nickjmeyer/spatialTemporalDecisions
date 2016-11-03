#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include "settings.hpp"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

Settings njm::sett;


const int Settings::numVals=2;

Settings::Settings(){
    isSetup = false;
    cleaned = 0;
    tick = std::time(NULL);
}

Settings::~Settings(){
    if(isSetup) {
        timeElapsed();
    }

    // compress data directory
    // if(!cleaned){
    //   int junk;
    //   junk=system(("cd " + srcDir + "; tar -cjf " + date + ".tar.bz2 "
    // 	 + date +"/").c_str());
    //   if(junk);// dummy statement to get -Wall off my back
    // }
}


void Settings::timeElapsed(){
    tock = std::time(NULL);

    std::stringstream timeInfo;
    timeInfo.str("");
    timeInfo.clear();
    timeInfo << std::endl << std::endl
             << "time elapsed: "
             << ((double)(tock-tick))/(3600.0)
             << " hours"
             << std::endl;
    std::cout << timeInfo.str();

    if(!cleaned){
        njm::toFile("\n\n\n",datDir + "/README.org");
        njm::toFile("* Time Elapsed",datDir + "/README.org");
        njm::toFile(timeInfo.str(),datDir + "/README.org");
    }
}


void Settings::setup(const std::string fileName, const std::string srcDir) {
    this->isSetup = true;
    this->fileName = fileName;
    this->srcDir = srcDir;


    timeStamp();

    // check input arguments
    // if(numInitVals!=numVals && numInitVals!=(numVals+1)){
    //   std::cout << "In Settings::Settings, "
    //       << "not the correct number of initial values.\n"
    //       << "Program requires " << numVals << " values.  "
    //       << "User provided " << numInitVals << " values.\n"
    //       << "Did you forget to count the executable as an argument?\n"
    //       << "Terminating program...\n";
    //   clean();
    //   exit(1);
    // }


    int arg=0;
    std::stringstream info;
    // fileName = initVals[arg++];

    char hostName[128];
    hostName[127] = '\0';
    gethostname(hostName,127);
    info << "hostName: " << hostName << "\n";

    info << "fileName: " << fileName << "\n";

    // srcDir=initVals[arg++];
    info << "srcDir: " << srcDir << "\n";

    datDir = srcDir + "/" + datDir;

    info << "datDir: " << datDir << "\n";
    info << "git-describe: " << buf.ptr << "\n";

    boost::filesystem::create_directories(datDir);

    // ask permission to proceed with program
    // std::string check;
    std::cout << info.str() << std::endl;
    // std::cout << "Are these settings correct? (y,n): ";

    // if(numInitVals==numVals) // did not pre specify the check?
    //   std::cin >> check;
    // else{
    //   check.assign(initVals[numVals]);
    //   std::cout << check << "\n";
    // }

    // if(check.compare("y")==0 || check.compare("Y")==0){
    std::cout << "Check complete.  "
              << "Proceeding with program...\n\n";
    njm::toFile("#+title: " + fileName,datDir + "/README.org");
    njm::toFile("#+author: Nick Meyer",datDir + "/README.org");
    njm::toFile("#+date: " + date,datDir + "/README.org");
    njm::toFile("\n\n\n",datDir + "/README.org");
    njm::toFile("* Run-time Information",datDir + "/README.org");
    njm::toFile(info.str(),datDir + "/README.org");
    //   // junk=system(("cp " + fileName + ".tar.bz2 " + datDir).c_str());
    //   // if(junk);// dummy statement to get -Wall off my back
    // }
    // else{
    //   std::cout << "Settings were not declared as correct.  "
    //       << "Terminating program...\n";
    //   clean();
    //   exit(0);
    // }

}


void Settings::clean(){

    // int junk;
    // std::string rm = "rm -rf " + datDir + " " + datDir + date + ".tar.bz2";
    // junk=system(rm.c_str());
    // if(junk);// dummy statement to get -Wall off my back
    boost::filesystem::remove_all(datDir);
    cleaned = 1;
}



void Settings::timeStamp(){
    time_t t=time(NULL);
    struct tm * lct;
    lct=localtime(&t);

    std::stringstream dirSS;

    dirSS << lct->tm_year+1900 << "-";

    if(lct->tm_mon+1 < 10)
        dirSS << "0";
    dirSS << lct->tm_mon+1 << "-";

    if(lct->tm_mday < 10)
        dirSS << "0";
    dirSS << lct->tm_mday << "-";

    if(lct->tm_hour < 10)
        dirSS << "0";
    dirSS << lct->tm_hour << "-";

    if(lct->tm_min < 10)
        dirSS << "0";
    dirSS << lct->tm_min << "-";

    if(lct->tm_sec < 10)
        dirSS << "0";
    dirSS << lct->tm_sec;

    date = dirSS.str();

    datDir="./" + dirSS.str();
    seconds=(int)t;
}



std::string Settings::datExt(int inDir, std::string beg, std::string end) const{
    std::stringstream ext;
    if(inDir)
        ext << datDir << "/";

    ext << beg;

    std::stringstream optSS;
    std::string opt;
    // begin options
    // end options

    opt=optSS.str();
    std::replace(opt.begin(),opt.end(),'.','p');
    ext << opt;

    ext << end;
    return ext.str();
}


std::string Settings::datExt(std::string beg, std::string end) const{
    return datExt(1,beg,end);
}


std::string Settings::datExt(int inDir, std::string end) const{
    return datExt(inDir,"",end);
}

std::string Settings::datExt(std::string end) const{
    return datExt(1,"",end);
}


std::string Settings::datExt(int inDir) const{
    return datExt(inDir,"","");
}


std::string Settings::datExt() const{
    return datExt(1,"","");
}


std::string Settings::srcExt(const std::string file){
    return srcDir + "/" + file;
}
