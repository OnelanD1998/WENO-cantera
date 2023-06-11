
#include <string>
#include <sstream>

#ifndef TIMETERM_H
#define TIMETERM_H

class timeTerm
{
private:
    double time_;
    double timeDelta_;
    dictionary timeDic_;
    double CFL_;
    int    overStep_;
    int     saveStep_;

public:
    timeTerm(/* args */);
    timeTerm(double startTime);
    ~timeTerm();

    dictionary& timeDic(){return timeDic_;}
    double& timeValue(){return time_;}
    double& timeDelta(){return timeDelta_;}
    double CFL(){return CFL_;}
    int overStep(){return overStep_;}
    int saveStep(){return saveStep_;}

    std::string timeName();
    void setTime(double timeNew){ time_ = timeNew; }
    void setTimeDelta(double timeNewDelta){ timeDelta_ = timeNewDelta; }
};

//--------------------------- constructor function ------------------//
timeTerm::timeTerm(/* args */)
:
time_(0),
timeDic_("timeDic")
{
    CFL_ = timeDic_.lookUp<double>("CFL");
    overStep_ = timeDic_.lookUp<double>("overStep");
    saveStep_ = timeDic_.lookUp<double>("saveStep");
}

timeTerm::timeTerm(double startTime)
:
time_(startTime),
timeDic_("timeDic")
{}

timeTerm::~timeTerm()
{
}


std::string timeTerm::timeName()
{
    using namespace std;
    stringstream conv;
    string  reci;
    conv << time_;
    conv >> reci;
    return reci;
}


#endif