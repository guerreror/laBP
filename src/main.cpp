/* 
 *  InvCoal
 *  Version:  labp_v17

 
 The core of this program simulates an Ancestral Recombination Graph (ARG) and constructs gene trees from it. Most of that ARG core simulation was written by Shane Pope and Mark Kirkpatrick (VIII-2008)
 
 This version, which simulates populations with an inversion polymorphism, was started by Rafael Guerrero VII-2009
 
 The results of the simulation (coalescent times at one site, for a sample of size 2) were confirmed using analytical models. The analytical work and simulation results are published in Guerrero, Rousset, & Kirkpatrick (Phil Trans Roy Soc B 2012). Additionally, the LD patterns between two sites have also been checked with analytical work (published in Rousset, Kirkpatrick & Guerrero, TPB 2014) and other simulations (Peischl, Koch, Guerrero & Kirkpatrick, Heredity 2013)
 
 Version notes:
 
 labp_v15: changed some of the hardcoded I/O to make the simulator easier to use in a pipeline. I'm also updating some of the code to be C++11 (barely). specifically, I will change the links to boost libraries that are now standard, and try to make the code cleaner.
 
 labp_v17: 
    -Added simple speciation (all carriers go to population 0, no forced coalescence, ancestral popSize is set to pop0)
    -"demography" events multiply all pop sizes by a coefficient
    -the inversion still has to be older than other events (speciation, demography)
    -Changed migration to allow for different population sizes
    -Added code to calculate Dxy, Fst
 
 */


// Includes from STL:
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::ostream;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include <string>
using std::string;
#include <random>
#include <chrono>

#include <memory>
using std::shared_ptr;
using std::weak_ptr;
using std::unique_ptr;


// Includes for our files:
#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "chromosome.h"
#include "sitenode.h"
#include "snptree.h"
#include "tajima.h"


// Global declarations:
std::random_device rd;
std::mt19937_64 gen(rd());		// Random generator declared globally. can be seeded in main.cpp and then used elsewhere


vector < vector <double> > buildMigMatrix(Parameters& p){
    // ======================================================================================
    // a simple stepping-stone migration matrix
    // Note that in other parts of the simulation migration assumes two populations only. Check migrateEvent()
    
    unsigned int nPops=static_cast<int>(p.paramData->popSizeVec.size());
    vector < vector <double> > mig_prob;
    mig_prob.resize(nPops);
    double m= p.paramData->migRate.at(0)/(2*p.paramData->totalPopSize);
    
    for(int i = 0; i < nPops; i++) mig_prob.at(i).resize(nPops);
    for(int i = 0; i < nPops; i++){
        for(int j = 0; j < nPops; j++){
            if (nPops==1) mig_prob.at(i).at(j)=1;
            else if (i==j) mig_prob.at(i).at(j)=1 - m;					// probability of staying in current pop
            else if (i==j-1) {
                mig_prob.at(i).at(j)= m/2;								// probability of migrating to neighboring pop at left
                if (i==0)  mig_prob.at(i).at(j)= m;						// probability of migrating if pop is at beginning of line
            }
            else if (i==j+1) { mig_prob.at(i).at(j)=m/2;				// probability of migrating to neighboring pop at left
                if (i==nPops-1)   mig_prob.at(i).at(j)= m;				// probability of migrating if pop is at end of line
            }
            else  mig_prob.at(i).at(j)=0;								// probability of migrating to any other pop
        }
    }
    
    /*
     // a set of N maps. In each map there are pairs of the following values:
     // <int> is the population
     // <float> is the cumulative probability of migration, to be compared against a random uniform number [0,1]
     //
     // Example: Set of four populations:: An individual in population 2 can migrate to pops {0,1,3} with prob {0,0.05, 0.05}
     // so the map for population 2 will look like this: { {0, 0}, {1, 0.95}, {2, 0.9}, {3, 1} }
     // For the use of these probabilities in migration, check out comments on migrate() function
     //
     //
     //
     */
    
    vector < map< double, int> > mig;
    mig.resize(nPops);
    for(int i = 0; i < nPops; i++){
        map< double, int> migrate;
        migrate[mig_prob.at(i).at(i)]= i ;
        double cumulative = mig_prob.at(i).at(i);
        for(int j = 0; j < nPops; j++){
            if (i!=j && mig_prob.at(i).at(j)!=0) {
                cumulative += mig_prob.at(i).at(j);
                migrate[cumulative]= j;
            }
        }
        mig.at(i)= migrate;
    }
    return mig_prob;
}

int main (int argc, const char * argv[]) {
    std::cerr<<"Seed: "<<rd()<<'\n';
    //// A static seed, useful for debuggin:
    //gen.seed(42U);
    
    // Timer start
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    
    unsigned int totalEvents=0;
    
    stringstream infile;
    stringstream ms_ss;
    stringstream sstat_ss;
    if(argc > 1 ) {
        infile <<"inLABP_"<< argv[1]<<".pars";
        ms_ss <<"outLABP_"<< argv[1]<<".sites";
        sstat_ss<<"outLABP_"<< argv[1]<<".stats";
    }
    else{
        infile <<"inLABP.pars";
        ms_ss <<"outLABP_" << rd()<<".sites";
        sstat_ss<<"outLABP_"<< rd()<<".stats";
    }
  
    Parameters params(infile.str().c_str());
    
    unsigned int nRuns=params.paramData->nRuns;
    unsigned nSites = params.paramData->n_SNPs;
 
    vector< vector<double> > mig_prob = buildMigMatrix(params);
    vector<double> outTime(params.paramData->n_SNPs, 0);

    std::ofstream msout;
    if(params.paramData->msOutput) {
        msout.open(ms_ss.str().c_str());
    }
    
    std::ofstream stout;
    if(params.paramData->msOutput) {
        stout.open(sstat_ss.str().c_str());
    }
    stout<<"ET totL S piTotal pi1 pi2 fst dxy tajD\n";

    double LDsum=0;
    int ticker=(int)(nRuns/10);
    
    if(ticker>0)std::cerr<<"Progress ticker (1 tick = 10% of runs): ";
    for (int timer=0; timer<nRuns; ++timer){
        if(ticker>0 && timer % ticker  == 0) std::cerr <<"+";
        
        params.setPhi();
        params.setSNPs();
        params.setCarriers();
        unsigned nCarriers = params.paramData->initChr.size();

        World * world = new World(params.getpData());
        
        // Run simulation until only a single carrier remains:
        while(!world->simulationFinished()){
            totalEvents += world->simulateGeneration(mig_prob);
        }
        
        
        vector< shared_ptr < ARGNode > > allNodes = world->getARGVec();
        vector<double>tempLD (params.paramData->n_SNPs,0);
        
        vector< SNPtree > trees;
        unsigned novar=0;
        vector< double > varPos;
        double lengthLastSite = 0;
        
        for(int k=0; k< nSites;++k){
            unsigned pos = k;

            if(!params.paramData->fixedS){pos = 0;} // if not doing fixed segregating sites, all sites are bases in a region with the same gene tree
            
            SiteNode geneTree = SiteNode(params.paramData->neut_site[pos], allNodes.at(allNodes.size()-1));
            SNPtree tmp(geneTree, nCarriers, params.paramData->theta);
            unsigned mutantCount = accumulate(tmp.SNPvalues.begin(), tmp.SNPvalues.end(), 0);
            
            if(mutantCount > 0){// this is a variable site, add to summary statistics
                trees.push_back(tmp);
                varPos.push_back(k);
            }
            else{ //No variation at this site. Keeping track
                ++novar;
            }
            lengthLastSite = tmp.totalLength/(double)params.paramData->totalPopSize;
            double totalmrca = geneTree.getTime()/(double)params.paramData->totalPopSize;
            outTime.at(k)+=  totalmrca;
            tempLD[k] = totalmrca;
        }
        
        double pi =  0;
        double piP1 =  0;
        double piP2 =  0;
        double dxy = 0;
        double fst = 0;
        
        double varSites = nSites - novar;

        for (int i=0; i< varSites; ++i){ //for every segregating site
            vector<unsigned> snps = trees.at(i).SNPvalues; //get sample at position
            vector< vector<unsigned> > split = split_by_population(snps, params.getSamplePerPop());
            double totalPi = heterozygosity(snps);
            pi+= totalPi;
            vector<double> pipop = pi_by_pop(split);
            piP1 += pipop[0];
            piP2 += pipop[1];
            fst += fst_nei(totalPi, pipop);
            dxy += calcdxy(split);
        }
        pi = pi /nSites;
        piP1 = piP1 /nSites;
        piP2 = piP2 /nSites;
        fst = fst /varSites;
        dxy = dxy / nSites;
        
        
        vector <double> posout =params.paramData->neut_site;
        if(!params.paramData->fixedS) posout = varPos;
        
        msout<<"// \nsegsites: "<<varSites<<"\npositions: ";
        for(int k=0; k< varSites;++k){ msout<<" "<<posout.at(k)<<" ";}msout<<'\n';
        for(int i=0; i<nCarriers;++i){for(int k=0; k < varSites;++k){msout<< trees[k].SNPvalues[i];} msout<<'\n';}
        
        double tempProd=tempLD[0]*tempLD[nSites-1];
        LDsum+=tempProd;
        
        stout << tempLD[0]<<" "<<lengthLastSite<<" "<<varSites <<" "<<pi<<" "<<piP1<<" "<<piP2<<" "<<fst<<" "<<dxy<<" "<<tajd(nCarriers, varSites, pi)<<'\n';

        // Clean up:
        //
        delete world;
    }
    
    std::cerr.precision(6);
    std::cerr<< "\nMean LD (E[T1,n]) = "<< LDsum/(double)nRuns <<'\n';

    if(!params.paramData->fixedS){ std::cerr <<"Mean TMRCA = "<< outTime.at(0)/nRuns <<'\n';}
    else{
        std::cerr <<"Mean TMRCA per site (E[T1]...E[Tn])\n";
        for(int k=0; k< nSites;++k){
            std::cerr<<" "<<outTime.at(k)/nRuns;
        }
        std::cerr<<"\n";
    }
    
    msout.close();
    stout.close();
    
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    
    return 0;
}
