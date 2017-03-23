/**
 * chromosome.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *	Started by Shane Pope, VIII-08;  modified by Mark Kirkpatrick, IX-08

 *  Modifications to simulate inversions by Rafael Guerrero 2008-2010
 *	This class is used to represent chromosomes that carry segments that are ancestral to
 *	those in the sample at time 0.  The data in Chromosome are
 *		• the context (pop, inversion allele)
 *		• a vector of chromosome segments
 *		• a vector of inverted chromosome segments
 *		• the descendant node of the chromosome
 *
 *	Functionality includes operations for recombination, migration, and coalescence.
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
	#include <algorithm>
		using std::max;
		using std::min;

   #include <memory>
		using std::shared_ptr;
		using std::weak_ptr;

// Includes from our files:
	#include "typedefs.h"
	#include "chromosome.h"
	#include "argnode.h"
	#include "ran_mk.h"




Chromosome::Chromosome(Context ctxt, vector<Segment> chromSegments,shared_ptr < ARGNode > desc ){
//
// Constructs a chromosome
//
	DBG("Chromosome:: Constructing Chromosome()...")	
	chromData.reset(new ChromosomeData);				// Makes chromData a pointer to the private data
	// Copy the context data:
	chromData->context = ctxt;
	// Copy the segments:
	chromData->chromosomeSegments.reserve( chromSegments.size() );
	copy(chromSegments.begin(), chromSegments.end(), 
			back_inserter(chromData->chromosomeSegments) );
	// Set the descendant ARG node:
	chromData->descendant = desc;
}

Chromosome::Chromosome(Context ctxt, vector<Segment> chromSegments, Segment invRange, shared_ptr < ARGNode >  desc ){
	//
	// Constructs a chromosome that carries an inversion
	//
	//
	DBG("Chromosome:: Constructing Chromosome() with inversion...")	
	
	chromData.reset(new ChromosomeData);				// Makes chromData a pointer to the private data
	// Copy the context data:
	chromData->context = ctxt;
	// This part makes sure that the
	//
	if(ctxt.inversion==1){
		DBG("segs before")
		for(int i=0; i<chromSegments.size();++i) {DBG( chromSegments.at(i).L)}
		vector<Segment> SegsRight = splitSegments(chromSegments, invRange.R, false);
		vector<Segment> SegsResidue = splitSegments(chromSegments, invRange.L, false);
		for(int i=0; i<SegsRight.size(); ++i) chromSegments.push_back(SegsRight[i]);
		
		DBG("segs after")
		for(int i=0; i<chromSegments.size();++i) {DBG( chromSegments.at(i).L)}
		DBG("in residue")
		for(int i=0; i<SegsResidue.size();++i) {DBG( SegsResidue.at(i).L)}
		
		chromData->chromosomeSegments.reserve( chromSegments.size() );
		chromData->invertedSegments.reserve( SegsResidue.size() );
		copy(chromSegments.begin(), chromSegments.end(), 
			 back_inserter(chromData->chromosomeSegments) );
		copy(SegsResidue.begin(), SegsResidue.end(), 
			 back_inserter(chromData->invertedSegments) );
		DBG("	chrom is inverted. Size of invSegs ("<<chromData->invertedSegments.size()<<") and Segs("<<chromData->chromosomeSegments.size()<<")\n")
	}
	
	else{
		chromData->chromosomeSegments.reserve( chromSegments.size() );
		copy(chromSegments.begin(), chromSegments.end(), 
			 back_inserter(chromData->chromosomeSegments) );
	}
    if(ctxt.inversion==0 && chromData->invertedSegments.size()!=0){std::cerr<<"Warning: construction of inverted chromosome in standard context)\n";}

	// Set the descendant ARG node:
	chromData->descendant = desc;
}

Chromosome::~Chromosome(){
//
// Destructor
//
	DBG("Chromosome:: Destructing...")

}



Chromosome::Chromosome( Chromosome& chr ){
//
// Copy constructor, needed to put in std::vector
//
	DBG("Chromosome:: Copy constructor called")	
	chromData.reset(new ChromosomeData);	
	// Copy the context data:
	chromData->context = chr.chromData->context;
	// Copy the segments:
	chromData->chromosomeSegments.reserve( chr.chromData->chromosomeSegments.size() );
	copy(chr.chromData->chromosomeSegments.begin(),
		 chr.chromData->chromosomeSegments.end(),
		 back_inserter(chromData->chromosomeSegments));
	
	if(chr.getInv()==1){
		chromData->invertedSegments.reserve( chr.chromData->invertedSegments.size() );
		copy(chr.chromData->invertedSegments.begin(),
			 chr.chromData->invertedSegments.end(),
			 back_inserter(chromData->invertedSegments));
	}
		
	// Set the descendant ARG node:
	chromData->descendant = chr.getDescendant();
}
	
void Chromosome::operator = ( Chromosome& chr ){
//
// Equals operator;  same result as copy constructor
//
	DBG("Chromosome:: = Operator called")	
	chromData.reset(new ChromosomeData);
	//Copy the context:
	chromData->context = chr.chromData->context;
	//Copy the segments:
	chromData->chromosomeSegments.reserve( chr.chromData->chromosomeSegments.size() );
	copy(chr.chromData->chromosomeSegments.begin(),
		 chr.chromData->chromosomeSegments.end(),
		 back_inserter(chromData->chromosomeSegments));	
	
	if(chr.chromData->context.inversion==1){
		chromData->invertedSegments.reserve( chr.chromData->invertedSegments.size() );
		copy(chr.chromData->invertedSegments.begin(),
			 chr.chromData->invertedSegments.end(),
			 back_inserter(chromData->invertedSegments));
	}
	
	// Set the descendant ARG node:
	chromData->descendant = chr.getDescendant();
}

shared_ptr<Chromosome> Chromosome::copy_values(Segment& inv_range)	{
    vector<Segment> segs=get_All_SegVec();
    assert(segs.size()!=0);
    shared_ptr<Chromosome> newChromosome(new Chromosome(getContext(), segs, inv_range, getDescendant()));
    
    return newChromosome;
}

void Chromosome::setInvSegs(vector<Segment> segs){
	chromData->invertedSegments.clear();	
	chromData->invertedSegments.reserve( segs.size() );
	copy(segs.begin(),
		 segs.end(),
		 back_inserter(chromData->invertedSegments));
}

void Chromosome::setStdSegs(vector<Segment> segs){
	chromData->chromosomeSegments.clear();
	chromData->chromosomeSegments.reserve( segs.size() );
	copy(segs.begin(),
		 segs.end(),
		 back_inserter(chromData->chromosomeSegments));
}

void Chromosome::clearInvSegs(){
	chromData->invertedSegments.clear();	
}

void Chromosome::setPopulation(unsigned short pop) {
	chromData->context.pop = pop;
}

void Chromosome::setInversion(unsigned short inv) {
	chromData->context.inversion = inv;
}

void Chromosome::setContext(Context c) {
	chromData->context= c;
}

void Chromosome::setDescendant(shared_ptr < ARGNode > desc) {
	//
	// Returns a pointer to the descendant ARG node
	//
	chromData->descendant = desc;
}

Context Chromosome::getContext() {
	return chromData->context;
}

unsigned short Chromosome::getPopulation() {
	return chromData->context.pop;	
}

unsigned short Chromosome::getInv() {
	return chromData->context.inversion;
}

shared_ptr < ARGNode > Chromosome::getDescendant() {
	return chromData->descendant;
}


vector<Segment> Chromosome::get_All_SegVec(){
	//
	// Returns vector of chromosome segments
	//
	// UPDATED JUNE2011, fixed bug
	
	vector<Segment> stdSegs=chromData->chromosomeSegments;
	vector<Segment> invSegs=chromData->invertedSegments;	
	vector<Segment> allSegs;
	
	DBG("get_allSegVec called. invSeg size= "<<invSegs.size()<<", StdSegs Size = "<<stdSegs.size())
	if(getInv()==0) return stdSegs;
	
	else {
		int i=0;
		if (invSegs.size()!=0 && stdSegs.size()!=0){
			while (stdSegs[i].R < invSegs[0].L && i<stdSegs.size()) {
				allSegs.push_back(stdSegs[i]); 
				++i;
				}
		}
		for(int j=0; j< invSegs.size(); ++j) allSegs.push_back(invSegs[j]);
		
		for(int k=i; k< stdSegs.size(); ++k) allSegs.push_back(stdSegs[k]);
		
		return allSegs;
	}
}

vector<Segment> Chromosome::get_Inv_Segs(){
	//
	// Returns vector of inverted segments
	//
	vector<Segment> segs;
	
	
	if(getInv()==0) return segs;
	
	if(getInv()==1) return chromData->invertedSegments;
	
    std::cerr<<"Unexpected end of function Chromosome::get_Inv_Segs()\n";
    return(segs);
}

vector<Segment> Chromosome::get_Std_Segs(){
	//
	// Returns vector of standard segments
	//
 return chromData->chromosomeSegments;
	
}

vector<Segment> Chromosome::translateSegs(const Segment& invRange){
	vector<Segment> invSegs=chromData->invertedSegments;
	vector<Segment> transSegs;
	
	for(int i=0; i<invSegs.size(); ++i){
		Segment s = invSegs[invSegs.size()-1-i];
		Segment t = { tr_pt(s.R, invRange), tr_pt(s.L, invRange)};
		transSegs.push_back(t);
	}
	return transSegs;
}

double tr_pt(double p, const Segment& invRange){
	return invRange.L + (invRange.R - p);
}

Segment Chromosome::getHomoRange(Segment& invRange) {
//
// Get the total length of the chromosome available for homokaryotypic recombination
//
	Segment range;
	unsigned long nSegs = chromData->chromosomeSegments.size();
	unsigned long nInvSegs = chromData->invertedSegments.size();
	
	
	if(nSegs == 0 && nInvSegs==0) {
		std::cerr << "Error in getHomoRange: the carrier is empty\n";
	}
	
	if(getInv()==1 && nInvSegs!=0){
										  
			vector<Segment> transSeg= translateSegs(invRange);
		
			if(nSegs!=0){
				
				range.L = min(chromData->chromosomeSegments[0].L, transSeg[0].L);
				range.R = max(chromData->chromosomeSegments[nSegs-1].R, transSeg[nInvSegs-1].R);
		
			}
			else{
				range.L = transSeg[0].L;
				range.R = transSeg[nInvSegs-1].R;			
			}
		}
	else{
		range.L =chromData->chromosomeSegments[0].L;
		range.R = chromData->chromosomeSegments[nSegs-1].R;			
	}	
	
	return range;

	}

double Chromosome::getHomoLength(Segment& invRange) {
	if(isEmpty()) return 0;
	Segment range= getHomoRange(invRange);
	return range.R-range.L;
}
	
double Chromosome::getHeteroLength(Segment& invRange) {
	//
	// Get the length of the chromosome available for heterokaryotipic recombination OUTSIDE OF THE INVERSION
	//
	//UPDATED Nov/10, RG. Fixed bug in length 
	
	unsigned int nSegs = chromData->chromosomeSegments.size();
	if(nSegs == 0) return 0;
	
	double l_length = invRange.L - chromData->chromosomeSegments[0].L; if( l_length<0) l_length=0;
	double r_length = chromData->chromosomeSegments[nSegs-1].R - invRange.R; if( r_length<0) r_length=0;
	
	return r_length + l_length;
}

double Chromosome::calcBreakpoint(Segment& invRange, bool hetero){
	// gives a breakpoint for recombination, allowing for heterokaryotypic events OUTSIDE OF THE INVERSION
	// UPDATED NOV/10 RG fixed bug in heterolength
	double breakpoint;

	if (hetero){
		
		if(chromData->chromosomeSegments[0].L<invRange.L) {	
			breakpoint = randreal(0, getHeteroLength(invRange)) + chromData->chromosomeSegments[0].L;
			if(breakpoint > invRange.L) breakpoint = invRange.R + (breakpoint-invRange.L) ;
		}
		else breakpoint = randreal(0, getHeteroLength(invRange)) + invRange.R;
			
	}
	else breakpoint = randreal(getHomoRange(invRange).L, getHomoRange(invRange).R);		// Compute the breakpoint
	return breakpoint;
}


bool Chromosome::isEmpty(){
	
	unsigned int nSegs = chromData->chromosomeSegments.size();
	unsigned int nInvSegs = chromData->invertedSegments.size();
	if(nSegs == 0 && nInvSegs==0) return true;
	return false;
}

bool Chromosome::carriesSite(double x){
    vector<Segment>segs=get_All_SegVec();
    for (int i=0; i< segs.size(); ++i){
        if(x >= segs[i].L && x <=segs[i].R) return true;
    }
    return false;
}

