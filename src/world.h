/* world.h
 *  InvCoal
 *  Version:  labp_v17
 *  Edited by Rafael Guerrero 2008-2016
 *
 * Data and methods to simulate evolution backwards.
 *
 *
 *
 */


#ifndef WORLD_
#define WORLD_


#include <memory>
    using std::shared_ptr;
	using std::weak_ptr;
    using std::unique_ptr;

// Includes from STL:
	#include <vector>
		using std::vector;
	#include <map>
		using std::map;
		
// Includes from our files:
	#include "typedefs.h"
	#include "argnode.h"
	#include "parameters.h"


// Forward declaration:
	class Chromosome;
	
	
// ==================================================================



class World {
 private:
	struct WorldData;									// This structure contains all the private data
	unique_ptr<WorldData> worldData;					// Pointer to private data
	cluster_t cluster;
	
public:
	World(shared_ptr<Parameters::ParameterData> p);		
	~World();
	void makecluster (unsigned long nPops);
	bool simulationFinished();				// Returns true if the simulation is finishe
    bool sitesCoalesced();
    bool timeCheck(double waiting, double Rate);
	double nGenerations();
	void Generation_pp();
	unsigned long totalNCarriers();					// Returns total number of carriers in the world
	unsigned long ctxtNcarriers(const Context ctxt); //	total number of carriers in a given context
	unsigned long popNcarriers(unsigned long pop);				//	total number of carriers in a given population
	unsigned short  migrateEvent(vector < vector< double> >& mig_prob, vector<double>& rate, double total);
	unsigned short  coalesceEvent(vector<double>& rate, double total);
	unsigned short  recombineEvent(vector<double>& rate, double total, bool hetero, bool gflux);
	unsigned short  doubrecEvent(vector<double>& rate, double total);
	unsigned short  simulateGeneration(vector< vector< double > >& mig_prob);
	unsigned short migratePoisson(vector < vector< double> >& mig_prob);
	unsigned short  coalescePoisson();
	unsigned short  recombine_all();
	vector< shared_ptr < ARGNode > >& getARGVec();
	double getFreqI(unsigned long pop);
	vector<double> getFreqbyPop(Context c);
    vector<double> getSizebyPop(Context c);
	Context getRecombContext(Context c, bool hetero);
	Context randCtx(unsigned long pop);
	void updateContexts(Context originCtx);
	void updateFreqs(vector<double> newF);
    void freqStepToLoss();
    void driftFreqStep();
    void forceAllCoal();
    void updateToNextEpoch();
    void demoChange();
    void speciation();
	shared_ptr< Chromosome> recomb_Wrap(shared_ptr<Chromosome> chrom, bool hetero, bool gflux);
    vector< shared_ptr < ARGNode > > initialARGnodes; //public vector of pointers to the original ARG nodes. For use in output of sample
    
	//CCD change imported from sexCoal
	//vector<double> getRecombBreakpoints();
	//vector< std::pair<double, shared_ptr < ARGNode > > >& getOutputSegments();

};

struct World::WorldData		// Structure with the private data for World, accessed by an opaque pointer
{
	double generation;								// Generation number, starting at 0
	unsigned int nPops;										// Number of populations in World
	unsigned int nClust;
	vector<int> nChromos;
    bool kingman;
    bool drift;
    vector<double> snpSites;
	Segment invRange;
	double phi;
	Context originCtx;
	vector<double> freq;
    double ancesterFreq;
	vector< unsigned int > popSize;
	vector< shared_ptr < ARGNode > > argNodeVec;					// Vector of pointers to all ARG nodes in the simulation
	shared_ptr< vector< vector < shared_ptr<Chromosome> > > > carriers;		// Pointers to arrays of pointers for carriers
    unsigned short current_epoch;
    vector< double > epoch_breaks;
    double epoch_Ncoef;
    bool epochs_over;
    vector<unsigned int> epochType;
    //CCD changes imported from sexCoal
    //vector<double> recomb_breakpoints;
    //vector< std::pair<double, shared_ptr < ARGNode > > > outputSegments;
    //double LeftTotalSegSize;
};



#endif /*WORLD_*/

