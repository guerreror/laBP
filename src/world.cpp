/**
 * World.cpp
 *  InvCoal
 *  Version:  labp_v17
 *  Edited by Rafael Guerrero 2008-2016
 *
 * This is the main class in our simulation
 * holds the vector of populations, the vectors of carriers per population, all parameter values
 * Methods for this class carry out the main events of the simulation (migration, recombination, coalescence)
 *
 *
 */

//#define DEBUGGER
// Includes from STL:
	#include <iostream>
		using std::cout;
		using std::endl;
	#include <vector>
		using std::vector;
	#include <map>
		using std::map;
	#include <math.h>
	
// Includes from our files:
	#include "chromosome.h"
	#include "ran_mk.h"
	#include "world.h"
	#include "argnode.h"

// ================================================================================================================


World::World(shared_ptr<Parameters::ParameterData> p){	
	// Constructor for World
	//

	DBG("Population:: Constructing Population(size, shared_ptr<Chromosome>)...")

	unsigned int nPops = p->popSizeVec.size();
	
	makecluster(nPops); // builds a map of all contexts
	
	unsigned int nClust=cluster.size();
	

	// Wipe clean WorldData:
	//
	worldData.reset(new WorldData);
	worldData->carriers.reset( new vector < vector< shared_ptr<Chromosome> > > );
	worldData->generation = 0;
	worldData->nPops = nPops;
	worldData->nClust = nClust;
    worldData->kingman=p->kingman;
	worldData->popSize=p->popSizeVec;
	worldData->drift=p->drift;
    worldData->snpSites = p->neut_site;
	worldData->invRange = p->invRange;
	worldData->phi=p->phi;
	worldData->freq=p->initialFreqs;

    worldData->freq.resize(nClust);
    for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter )
    {    worldData->freq[(*iter).second]= 1 - p->initialFreqs[(*iter).first.pop];
        if((*iter).first.inversion==1)
            worldData->freq[(*iter).second]= p->initialFreqs[(*iter).first.pop];
    }
    
    worldData->current_epoch = 0;
    
    if(p->demography.at(0)==1) {
        worldData->epoch_Ncoef = p->demography.at(2);
        if(p->speciation.at(0)==1){
            worldData->ancesterFreq = p->speciation.at(2);
            
            if(p->speciation.at(1)<p->demography.at(1)){
                worldData->epoch_breaks.push_back(p->speciation.at(1));
                worldData->epoch_breaks.push_back(p->demography.at(1));
                worldData->epochType.push_back(1);
                worldData->epochType.push_back(2);
            }
            else{
                worldData->epoch_breaks.push_back(p->demography.at(1));
                worldData->epoch_breaks.push_back(p->speciation.at(1));
                worldData->epochType.push_back(2);
                worldData->epochType.push_back(1);
            }
        }
        else{
            worldData->epoch_breaks.push_back(p->demography.at(1));
            worldData->epochType.push_back(2);
        }
    }
    else{
        worldData->epoch_Ncoef =1;
        if(p->speciation.at(0)==1){
            worldData->ancesterFreq = p->speciation.at(2);
            worldData->epoch_breaks.push_back(p->speciation.at(1));
            worldData->epochType.push_back(1);
            }
        }
    
    if(p->inv_age>0){
        worldData->epoch_breaks.push_back(p->inv_age);
        worldData->epochType.push_back(0);
    }
    
    worldData->epochs_over = false;
    if(worldData->epoch_breaks.size()==0) worldData->epochs_over = true;

    Context originCtx(0,0);
	worldData->originCtx=originCtx;
	// Copy the population sizes and resize vectors to the number of contexts:
	//
	worldData->carriers->resize(nClust);
	
	DBG("n of contexts "<<cluster.size())
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) 
	{
		DBG((*iter).first.inversion<<" inversion, "<<(*iter).first.pop<<" pop is: "<< (*iter).second<<" freq = "<<worldData->freq.at((*iter).second));
	}
	
	
	vector<int> initNCarriers;
	initNCarriers.resize(nClust);
	shared_ptr< ARGNode > nulle ;
	
    for(int i = 0; i < p->initChr.size(); ++i){
        assert(!p->initChr.at(i)->isEmpty());
        shared_ptr<Chromosome> chr= p->initChr.at(i)->copy_values(worldData->invRange);
        shared_ptr < ARGNode > newNode ;newNode.reset( new ARGNode(worldData->argNodeVec.size(), chr) );// Make a terminal ARG node
        chr->setDescendant(newNode);
        worldData->carriers->at(cluster[chr->getContext()]).push_back( chr );					// Adds this chromosome to the carriers array
        worldData->argNodeVec.push_back(newNode);							// Add the new node to the vector of ARG nodes
        
        initialARGnodes.push_back(newNode); //added X-13, saved for alternative output
    }
    
    
    worldData->generation=1.0;
} // End constructor for World



World::~World(){
	DBG("Population:: Destructing...")
}



bool World::simulationFinished(){
    //
    // Determine if only one carrier remains in the world
    //
    if(totalNCarriers()==1) return true;

//C.Cheng, change end condition
//	if( worldData->LeftTotalSegSize <= 0 || totalNCarriers()==1) return true;
//	if (worldData->recomb_breakpoints.size() == 0) return true;

    return false;
}

bool World::sitesCoalesced(){
    //
    // Determine if all sites of interest have coalesced
    //
    vector<int> carrierCounter(worldData->snpSites.size(),0);

	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
        for(int j=0; j< worldData->carriers->at((*iter).second).size(); ++j){
            for(int s=0; s< worldData->snpSites.size(); ++s){
                double site= worldData->snpSites[s];
                if(worldData->carriers->at((*iter).second)[j]->carriesSite(site)) ++carrierCounter[s];
            }
        }
    }
    bool out = true;
    for (int i=0; i< carrierCounter.size();++i){
        if(carrierCounter[i]>1) out= false;
    }
    return out;
}

void World::forceAllCoal(){
 
	vector<shared_ptr<Chromosome> > sons;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
			unsigned long n=ctxtNcarriers((*iter).first);
			
			for(int j=0; j<n;++j){
				shared_ptr<Chromosome> chr=worldData->carriers->at((*iter).second).at(j);
                chr->setStdSegs(chr->get_All_SegVec());
				chr->setContext(worldData->originCtx);
				sons.push_back(chr);
            }
			worldData->carriers->at((*iter).second).clear();
	}
    
	if(sons.size()>0){
		shared_ptr<Chromosome> cain=sons.at(0);
		
		shared_ptr < ARGNode > cainNode ;
		cainNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, worldData->generation) );
		
		for(int i=1; i<sons.size(); ++i){
			shared_ptr<Chromosome> abel=sons.at(i);
			shared_ptr < ARGNode > newNode ;
			newNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, abel, worldData->generation) );
			worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation
			cain->merge( abel, newNode); // merging includes setting decendant to new ARGNode.
		}
		
		worldData->carriers->at(cluster[worldData->originCtx]).push_back(cain);
	}

}


unsigned long World::ctxtNcarriers(const Context ctxt){
	unsigned long total = worldData->carriers->at(cluster[ctxt]).size();
	return total;
}


unsigned long World::popNcarriers(unsigned long pop){
// Returns the number of carriers in all contexts of population pop	
	unsigned long total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop)
			total = worldData->carriers->at((*iter).second).size() +total;
	}
	return total;
	
}

unsigned long World::totalNCarriers(){
//
// Returns the total number of carriers in the world
//
	unsigned long total = 0;
	for(unsigned long i = 0; i < worldData->nClust; i++)
	{	total = worldData->carriers->at(i).size() +total;
		}
	return total;
}



double World::nGenerations(){
	return worldData->generation;
}

void World::Generation_pp (){
	worldData->generation++;
}

void World::makecluster (unsigned long nPops){
	
	int iter=0;
	for (int i=0; i<nPops; ++i){			//loop through populations
		for (int j=0; j < 2; ++j){			//loop through inversion states (I,S)	//loop through alleles in siteA

				Context c(i, j); 
				cluster[c]=iter;
				++iter;
			
			}
		}

}

vector< shared_ptr < ARGNode > >& World::getARGVec(){
	return worldData->argNodeVec;
}

double World::getFreqI(unsigned long pop){
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (1 ==(*iter).first.inversion && pop==(*iter).first.pop)
			total+= worldData->freq.at((*iter).second);
	}
	return total;
	
}

vector<double> World::getFreqbyPop(Context c){
	vector<double> freqs (worldData->nPops,0);
	Context geno=c;
	for (int i=0; i<worldData->nPops;++i){
		geno.pop=i;
		freqs.at(i)=worldData->freq.at(cluster[geno]);
	}
	return freqs;
}

vector<double> World::getSizebyPop(Context c){
    vector<double> sizes (worldData->nPops,0);
    Context geno=c;
    for (int i=0; i<worldData->nPops;++i){
        geno.pop=i;
        sizes.at(i)=worldData->freq.at(cluster[geno]) * worldData->popSize.at(geno.pop);
    }
    return sizes;
}


Context World::getRecombContext(Context c, bool hetero){
	// this function gives the context of a chromosome that will recombine with a carrier
	// given that the recombination is homokaryotypic or heterokaryotipic
	//  the function randomly selects one of the possible genotypes within standard or inverted contexts.
	// 
	
	//if(hetero)cout<<hetero<<'\n';
	double f= getFreqI(c.pop);
	int inv=1;
	if ((c.inversion==1 && hetero) || (c.inversion==0 && !hetero)){ 
		f=1-f;
		inv=0;
	}
	double draw = randreal(0, f);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (inv ==(*iter).first.inversion && c.pop==(*iter).first.pop){
			total+= worldData->freq.at((*iter).second);
			if(draw<total) {
					return (*iter).first;
			}
		}
	}
	cout<<"Error in getRecombContext. hetero= "<<hetero<<", f= "<<f<<", inv= "<<inv<<", total= "<<total<<"\n";
	return c;
}

void World::updateContexts(Context originCtx){
	// this function will transfer all carriers on inverted contexts to standard context. It's meant to be done 
	// at the origin of the inversion
	
	vector<shared_ptr<Chromosome> > sons;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (1 ==(*iter).first.inversion && ctxtNcarriers((*iter).first)>0 ){
			unsigned long n=ctxtNcarriers((*iter).first);
			
			for(int j=0; j<n;++j){
				shared_ptr<Chromosome> chr=worldData->carriers->at((*iter).second).at(j);
				chr->setStdSegs(chr->get_All_SegVec());
				chr->clearInvSegs();
				chr->setContext(originCtx);
				sons.push_back(chr);
			}
			worldData->carriers->at((*iter).second).clear();
		}
		
	}
	if(sons.size()>0){
		shared_ptr<Chromosome> cain=sons.at(0);
		
		//UPDATE NOV/10 RG included these two lines to make sure the node itself changes context when switched to Standard (even without coalescence)
		shared_ptr < ARGNode > cainNode ;
		cainNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, worldData->generation) );
		
		for(int i=1; i<sons.size(); ++i){
			shared_ptr<Chromosome> abel=sons.at(i);
			shared_ptr < ARGNode > newNode ;
			newNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, abel, worldData->generation) );
			worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation	
			cain->merge( abel, newNode); // merging includes setting decendant to new ARGNode.
		}
		
		worldData->carriers->at(cluster[originCtx]).push_back(cain);
	}
}

void World::updateFreqs(vector<double> newF){
	worldData->freq=newF;
}
	
Context World::randCtx(unsigned long pop){
	double draw = randreal(0, 1);
	Context c(0,0);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop){
			total+= worldData->freq.at((*iter).second);
			if(draw<total) {
				return (*iter).first;
			}
		}
	}
	cout<<"Error in randCtx\n";
	return c;
	
}
