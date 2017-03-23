/**
 *  parameters.h
 *  InvCoal
 *  Version:  labp_v17
 *
*/

#ifndef PARAMETERS_H_
#define PARAMETERS_H_


#include <memory>
    using std::shared_ptr;
    using std::unique_ptr;

    using std::vector;
#include <numeric>
    using std::accumulate;

// Includes from our files:
	#include "typedefs.h"


// Forward declaration:
	class Chromosome;

class Parameters {
  
 public:
	// private:
	struct ParameterData;
	shared_ptr<ParameterData> paramData;
	
	Parameters();
	Parameters(const char* insstring);
	~Parameters();
	vector<unsigned int> getPopulationSizes();
    vector<unsigned int> getSamplePerPop();
	vector<shared_ptr<Chromosome> > getChromVec();
	shared_ptr<ParameterData> getpData();
    void setPhi();
    void setCarriers();
    void setSNPs();
    void setTimeStamp(double t);
	
};

struct Parameters::ParameterData{
	//
	// Private data for Parameters
	//
    unsigned int nRuns;
	vector<unsigned int> popSizeVec;	    	// Vector with population sizes (numbers of chromosomes) in each population
	vector<double> migRate;
	vector<double> initialFreqs; //Vector of frequencies of inversion (one per population)
    vector<double> speciation;
    vector<double> demography;
    unsigned int inv_age;
    vector<double> phi_range;
    Segment invRange;
    double phi;
    double theta;
    double BasesPerMorgan;
    bool kingman;
    bool drift;
    bool msOutput;
    unsigned int n_SNPs;
    bool randSNP;
    bool randPhi;
    bool fixedS;
    vector<double> snpRange;
    vector<double> snpPositions;
    vector<vector<int> > nCarriers;
	vector<shared_ptr<Chromosome> >initChr;		// vector of initial carriers
	vector<double> neut_site;
    double timestamp;
    unsigned long totalPopSize;
};


#endif /*PARAMETERS_H_*/


