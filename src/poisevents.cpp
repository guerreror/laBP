/*
 *  poisevents.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *  Created by Rafael Guerrero on 6/9/09.
 *  These are methods for the World class, used during simulation of a generation-by-generation coalescent process.
 *  
 *
 */
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <math.h>
using std::pair;
using std::make_pair;
#include <algorithm>
using std::count;

#include <boost/math/special_functions/binomial.hpp>
using boost::math::binomial_coefficient;

// Includes from our files:
#include "chromosome.h"
#include "ran_mk.h"
#include "world.h"
#include "argnode.h"



unsigned short World::migratePoisson(vector < vector< double> >& mig_prob){
	//UPDATE DEC10, RG: Migration rates adjusted to account for nonconservative migration
	//typical of scenarios like local adaptation
	//originally, the migration rate per individual was m
	// corrected, backward migration rate S1->S2 = m*q1/q2 where qi is the size of the context i
	//
	//The function was completely rewritten.
	
	// Initialize things:

    vector<pair<shared_ptr<Chromosome>,int> > from_to;
	vector<double> rate;
	
	
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {					
		
		unsigned int k=ctxtNcarriers(i->first);
		int pop= i->first.pop;																	
		int c= i->second;
		vector<double> mfreq = getSizebyPop(i->first);
		vector<double> migs;
		
		
		for(int l=0; l<worldData->nPops;++l){
			if(l!=pop){
				double b_mig = mig_prob.at(pop).at(l) * mfreq.at(l)/mfreq.at(pop); 
                unsigned int ev = randpois((double)k * b_mig);
				if(ev > k) ev = k;
				vector<int> carr_idx;
				for (int u=0; u< worldData->carriers->at(c).size(); ++u) carr_idx.push_back(u);
				random_shuffle(carr_idx.begin(), carr_idx.end());
				
				for(int o=0;o<ev;++o){
					shared_ptr<Chromosome> who = worldData->carriers->at(c).at(carr_idx.at(o));
					from_to.push_back( make_pair(who,l));
				}
			}
		}
		
	}
	
	
	for (int i=0; i<from_to.size(); ++i){
		
		shared_ptr<Chromosome> chrom = from_to.at(i).first;	// Get the migrant carrier
		
		int c = cluster[chrom->getContext()];
		worldData->carriers->at(c).erase(remove(worldData->carriers->at(c).begin(), worldData->carriers->at(c).end(), chrom), worldData->carriers->at(c).end());
		
		chrom->setPopulation(from_to.at(i).second);									// change the context of this carrier to its new population
		
		int newC= cluster[chrom->getContext()];	
		worldData->carriers->at(newC).push_back( chrom );								// add the carrier to the vector of its new population
		
		shared_ptr < ARGNode >  newNode ;
		newNode.reset( new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation ) );
		worldData->argNodeVec.push_back(newNode);						// Add a pointer to this node to the vector of all nodes in	the simulation
		chrom->setDescendant(newNode);
		
	}
	
	
	return from_to.size();
}	//********************************************* end of migratePoisson

unsigned short World::coalescePoisson (){
	// Written by Rafael Guerrero, VI-09
	//
	
	int totalEvents=0;
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {	
		unsigned int k=ctxtNcarriers(i->first);
		unsigned long clustID=i->second;
		int nLocalChromos = worldData->freq.at(clustID) * worldData->popSize.at(i->first.pop);
		
		if(k > 1){												
			vector<shared_ptr< Chromosome> > coalescers;
			coalescers.resize(0);
			
			for(int carrier = 0; carrier < k; carrier++){
				
				shared_ptr<Chromosome> who = worldData->carriers->at(clustID).at(carrier);
				
				int ranParent = randint(0, nLocalChromos-1);
				if( ranParent < coalescers.size()){			//coalescence happened
					totalEvents++;
					shared_ptr<Chromosome> receiver = coalescers.at(ranParent);
					
					// Make a new ARG node: 
					//
					shared_ptr < ARGNode >  newNode ;
					newNode.reset( new ARGNode(worldData->argNodeVec.size(), receiver, who, worldData->generation ) );
					worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation
					
					// Merge the coalescing chromosomes
					//
					receiver->merge( who, newNode );
				}
				else{
					coalescers.push_back(who);
				}
			}
			
			worldData->carriers->at(clustID) = coalescers;
		}
	}
	
	return totalEvents;
}//**********************************************End of CoalescePoisson


unsigned short World::recombine_all(){
	//
	//	goes through list of carriers and asks if recombination happened
	//
	//	Written by RG VI-09 // UPDATE III-11 by RG	//
    //  Modified into recombine_all() by RG (III-16), does gene flux, heterokaryotypic, and homokaryotipic events.
    
	int events=0;
	
	vector< shared_ptr < Chromosome> > recombinants;
	
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
		unsigned long clustID=i->second;
		
		for(int carrier =0 ; carrier < worldData->carriers->at(clustID).size(); ++carrier){
			shared_ptr< Chromosome > chrom= worldData->carriers->at(clustID).at(carrier);
            
            double  freqI = getFreqI(i->first.pop);
            if (i->first.inversion==0) freqI= 1- freqI;

			double gfRate = (1-freqI)* worldData->phi;
            double hetRate = (1- freqI)* chrom->getHeteroLength(worldData->invRange);
            double homRate = freqI * chrom->getHomoLength(worldData->invRange);
	
            double x = randreal(0, 1);
			
			if (x < (gfRate + hetRate + homRate)){
				events++;
				
                bool gflux= false;
                bool hetero = false;
                
                if(x < gfRate) {
                    gflux = true;
                }
                else if(x < gfRate+hetRate){
                    hetero=true;
                }
                
				// creat new chromosome
				shared_ptr<Chromosome> chrom2 = recomb_Wrap(chrom, hetero, gflux);
				
				CDBG("After recomb, chr has "<<chrom->get_All_SegVec().size()<<" segs, and chr2 "<< chrom2->get_All_SegVec().size())
				
				if(!chrom->isEmpty()) {
					recombinants.push_back(chrom);
					CDBG("		chr has segs")
					//cout<<cluster[chrom->getContext()]<<'\n';
				}
				
				if(!chrom2->isEmpty()) {
					recombinants.push_back(chrom2);	//	add chrom2 to the world only if it holds segments
					CDBG("		chr2 has segs")
					//cout<<cluster[chrom2->getContext()]<<'\n';
				}
				
			}
			else{
				recombinants.push_back(chrom);
				CDBG("no recombination. pushing chrom to next gen")
			}
		}
	}
	
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
		worldData->carriers->at(i->second).clear();
	}
	
	for(int i=0; i< recombinants.size(); ++i){
		shared_ptr< Chromosome > chrom= recombinants[i];
		Context c= chrom->getContext();
		worldData->carriers->at(cluster[c]).push_back(chrom);
		
		CDBG("		carrier added in "<<cluster[c]<<", carriers ="<<worldData->carriers->at(cluster[c]).size())
	}
	
	return events;
}
