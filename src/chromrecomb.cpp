/*
 *  chromrecomb.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *  Created by Rafael Guerrero on 10/3/10.
 *
 */
// Debugging defines/includes
//#define DEBUGGER
//#define RCDBGR
//#define NDEBUG  // this one turns OFF assert fn
//#define CDBGR

#include <assert.h>

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

vector<Segment> splitSegments(vector<Segment>& old, double breakpoint, bool inverted){
	//
	// Erases segments on this chromosome to the right of breakpoint, and returns the vector
	//	of segments that were erased.  If the breakpoint falls within a segment, the part of
	//	that segment to the right of the breakpoint is erase, and that part of the segment
	//	is included in the vector of segments that are returned.
	//	If inverted==true, then right and left products are switched.
	//
	//
	
	SRCDBG("splitting at "<<breakpoint )
	
	vector<Segment> newSegVec;										// New vector of segments that lie to right of breakpoint
	unsigned long nSegs = old.size();
	for(int i = 0; i < nSegs; ++i){									// Loop through segments carried by this chromosome
		
		if(breakpoint < old[i].R){		// Breakpoint falls to the left or within this segment
			SRCDBG("		"<<breakpoint <<" < "<<old[i].R)
			int within=0;
			if(old[i].L < breakpoint)
			{												// Breakpoint falls within this segment
				within=1;
				SRCDBG("			breakpoint within segment, "<<old[i].L <<" < "<<breakpoint)	
				
				// Add partial segment to the right of breakpoint to the new segment vector, and 
				//	erase it from this chromosome:
				//
				Segment newSeg = {breakpoint, old[i].R};
				SRCDBG("			New seg "<<newSeg.L<<" -- "<<newSeg.R)
				newSegVec.push_back(newSeg);
				old[i].R = breakpoint;
			}
			else {
				SRCDBG("			adding seg "<< i)
				newSegVec.push_back(old[i]);
			}
			
			for(int j = i+1; j < nSegs; j++)	
			{															// Add remaining segments to the new segment vector
				SRCDBG("			adding seg "<< j)
				newSegVec.push_back(old[j]);	
			}
			old.resize(i + within);										// Erase remaining segments from this chromosome
	
			if(inverted){
				SRCDBG("Inverted segments. Switching content.")
				vector<Segment> temp = old;
				old= newSegVec;
				newSegVec=temp;
			}
			SRCDBG("newSeg Size= "<<newSegVec.size()<< ", old size= "<<old.size())
			return newSegVec;											// Return
			
		} // end breakpoint falls to the left or within this segment
	} // end loop through segments
	
	newSegVec.resize(0);		// Breakpoint falls to the right of the rightmost segment carried by this chromosome
	SRCDBG("newSeg empty")
	
	
	if(inverted){
		SRCDBG("Inverted segments. Switching content.")
		vector<Segment> temp = old;
		old= newSegVec;
		newSegVec=temp;
	}
	SRCDBG("newSeg Size= "<<newSegVec.size()<< ", old size= "<<old.size())
	return newSegVec;
}


shared_ptr<Chromosome> Chromosome::recombine(double breakpoint, shared_ptr< ARGNode > descNode, Context newContext, Segment& invRange) {
	//
	// Recombine this chromosome at breakpoint.  The chromosome that is operated on keeps the
	//	segments to the left of the breakpoint. 
	//	and returs a new chromosome with the righthand segments as the value of the function.
	//
	//
	assert(!isEmpty() );
	if(getInv()==0) assert(chromData->invertedSegments.size()==0);
	
	bool hetero= getInv()!=newContext.inversion;
	bool within=breakpoint>invRange.L && breakpoint<invRange.R;
	bool two_inverted= getInv()==1 && newContext.inversion==1;

	assert(!(hetero&&within));
	
	vector<Segment> newSegs= splitSegments(chromData->chromosomeSegments, breakpoint, false);
	vector<Segment> invSegs;
	
	RCDBG("\n\nRecomb btwn inversion:"<<getInv()<<" and "<<newContext.inversion)
	RCDBG("break at "<<breakpoint)
	RCDBG("NewSegs "<<newSegs.size()<<" and this SegVec "<<chromData->chromosomeSegments.size())
	RCDBG("within = "<<within)
	
	if(breakpoint < invRange.L) {
		if(hetero)	{int t=getInv(); setInversion(newContext.inversion); newContext.inversion=t;}
		vector<Segment> temp=chromData->invertedSegments;chromData->invertedSegments=invSegs; invSegs=temp;
	}
		
	if(two_inverted && within){
		// the calculation of the breakpoint gives a 'translated' breakpoint. 
		//that is why, if the break is inside of the inversion we re-translate it
		RCDBG("two_inverted. Segs before: ")
		//for (int i=0; i<chromData->invertedSegments.size();++i) cout<<"		"<<chromData->invertedSegments.at(i).L;
		//cout<<'\n';
		
		breakpoint=tr_pt(breakpoint, invRange);
		RCDBG("new bp ="<< breakpoint)
		
		invSegs = splitSegments(chromData->invertedSegments, breakpoint, true);
		
		RCDBG("Segs after: ")
		//for (int i=0; i<chromData->invertedSegments.size();++i) cout<<"		"<<chromData->invertedSegments.at(i).L;
		//cout<<'\n';
		RCDBG("Segs in invSegs: ")
		//for (int i=0; i<invSegs.size();++i) cout<<"		"<<invSegs.at(i).L;
		//cout<<'\n';
	}
	// Make a new chromosome that carries the segments to the right of the breakpoint:
	//
	
	shared_ptr<Chromosome> newChromosome(new Chromosome(newContext, newSegs, invRange, descNode));		
	newChromosome->setInvSegs(invSegs);		

	if(!isEmpty() && !newChromosome->isEmpty()){	// if both carriers have material, update ARGNode
		//cout<<'\n';
		newChromosome->setDescendant(descNode);
		chromData->descendant = descNode;								// Update the descendant ARG node for this chromosome:
	}
	else{	// else, the old ARGNode is used and the recomb event will be "ignored". Empty chromosomes will be deleted in recombineEvent.
		newChromosome->setDescendant(chromData->descendant);
	}
	
	RCDBG("Resulting in :"<<getInv()<<getPopulation()<<" and "<<newContext.inversion<<newContext.pop)
	RCDBG("Empty in old ("<<(int)isEmpty()<<") or new("<<(int)newChromosome->isEmpty()<<")")
	

	return newChromosome;
}

shared_ptr<Chromosome> Chromosome::doubrecombine(double breakpoint1, double breakpoint2, shared_ptr< ARGNode > descNode, Context newContext, Segment& invRange) {
	//
	// Recombine this chromosome at breakpoint.  The chromosome that is operated on keeps the
	//	segments to the left of the breakpoint. 
	//	and returs a new chromosome with the righthand segments as the value of the function.
	//
	assert(!isEmpty());
	assert(breakpoint1<breakpoint2);
	
	vector <Segment> old;
	if(getInv()==1) old=chromData->invertedSegments;
	else old=chromData->chromosomeSegments;
	DBG("old Segs from inv="<< getInv())
	
	// Make a new segment vector that carries the segments to the right of the rightmost breakpoint, then cut again at the
	// left breakpoint, paste the left and right segments (leaving the middle as Residue):	
	vector<Segment> SegsRight = splitSegments(old, breakpoint2, false);
	vector<Segment> SegsResidue = splitSegments(old, breakpoint1, false);
	for(int i=0; i<SegsRight.size(); ++i) old.push_back(SegsRight[i]);
	//update this chrom's segments
	if(getInv()==1) chromData->invertedSegments=old;
	else chromData->chromosomeSegments=old;
	DBG("old Size= "<< old.size()<<" and new = "<<SegsResidue.size() )

	
	// Make a new chromosome that carries the segments between the breakpoints:
	//
	shared_ptr<Chromosome> newChromosome(new Chromosome(newContext,SegsResidue, invRange, descNode));


if(!isEmpty() && !newChromosome->isEmpty()){	// if both carriers have material, update ARGNode
	newChromosome->setDescendant(descNode);
	chromData->descendant = descNode;								// Update the descendant ARG node for this chromosome:
}
else{	// else, the old ARGNode is used and the recomb event will be "ignored". Empty chromosomes will be deleted in recombineEvent.
	newChromosome->setDescendant(chromData->descendant);
}
	
	DBG("Resulting in :"<<getInv()<<getPopulation()<<" and "<<newContext.inversion<<newContext.pop)
	DBG("Empty in old ("<<(int)isEmpty()<<") or new("<<(int)newChromosome->isEmpty()<<")")

return newChromosome;
}

void Chromosome::merge(shared_ptr <Chromosome> chrom, shared_ptr < ARGNode > descNode) {
	//
	// Coalesces this chromosome with chrom.  
	//	Sets the descendant to the new descNode
	//cout<<"Merge ";for(int i=0;i<get_Inv_Segs().size();++i)cout<<get_Inv_Segs().at(i).L<<" "; 
	//cout<<"and ";for(int i=0;i<chrom->get_Inv_Segs().size();++i)cout<<chrom->get_Inv_Segs().at(i).L<<" ";cout<<'\n';
	//cout<<"with Std ";for(int i=0;i<get_Std_Segs().size();++i)cout<<get_Std_Segs().at(i).L<<" "; 

	chromData->chromosomeSegments= combineDoubletVectors(chromData->chromosomeSegments, chrom->chromData->chromosomeSegments);
	if(getInv()==1)chromData->invertedSegments = combineDoubletVectors(chromData->invertedSegments, chrom->chromData->invertedSegments);
	
	//cout<<"Result ";for(int i=0;i<get_Inv_Segs().size();++i)cout<<get_Inv_Segs().at(i).L<<" "; 
	//cout<<" ----- ";for(int i=0;i<chrom->get_Inv_Segs().size();++i)cout<<chrom->get_Inv_Segs().at(i).L<<" ";cout<<'\n';
	//cout<<"with Std ";for(int i=0;i<get_Std_Segs().size();++i)cout<<get_Std_Segs().at(i).L<<" "; 

	chromData->descendant = descNode;
}


vector<Segment> Chromosome::combineDoubletVectors(vector<Segment> vector1, vector<Segment> vector2){
	//
	// Add a doublet to a doublet vector and fix overlapping etc.
	//
	vector<Segment> outputVector;
	
	if( vector1.size() == 0 && vector2.size() == 0 );	// Both vectors are empty (never happen hopefully)
	else if(vector1.size() == 0)						// vector1 is empty, copy vector2
		copy(vector2.begin(), vector2.end(), back_inserter(outputVector));
	else if(vector2.size() == 0)						// vector2 is empty, copy vector1
		copy(vector1.begin(), vector1.end(), back_inserter(outputVector));
	else{												// Neither vector is empty
		
		//This is a "state machine" with 7 states.
		//The first four states correspond to different kinds of comparisons between vector1 and vector2.
		// State 0: vector1's L compared to vector2's L
		// State 1: vector1's R compared to vector2's R
		// State 2: vector1's R compared to vector2's L
		// State 3: vector1's L compared to vector2's R 
		// State 4: vector1 has no more segments, copy over vector2
		// State 5: vector2 has no more segments, copy over vector1
		// State 6: both vectors empty, we are finished copying
		
		//Vector indices:
		unsigned long indexVec1 = 0;
		unsigned long indexVec2 = 0;
		unsigned long indexVecOut = 0;
		
		unsigned long STATE = 0;
		bool finished = false;
        
		while(!finished){
			switch (STATE){
				case 0:{ //vector1 -> L, vector2 -> L
					if(indexVec1 == vector1.size()){
						STATE = 4;
					}
					else if(indexVec2 == vector2.size()){
						STATE = 5;
					}
					else{
						Segment tempDoublet;
						outputVector.push_back(tempDoublet);
						indexVecOut = outputVector.size()-1;
						
						if(vector1.at(indexVec1).L == vector2.at(indexVec2).L){
							outputVector.at(indexVecOut).L = vector1.at(indexVec1).L;
							STATE = 1; //1 = R, 2 = L
		            	}
			            else if(vector1.at(indexVec1).L < vector2.at(indexVec2).L){
							outputVector.at(indexVecOut).L = vector1.at(indexVec1).L;
							STATE = 2; //1 = R, 2 = L
			            }
			            else{
			                outputVector.at(indexVecOut).L = vector2.at(indexVec2).L;
			                STATE = 3; //1 = L, 2 = R
			            }
					}
		            break;//case 0
	        	}
					
		        case 1:{	//1 = R, 2 = R
		            if(vector1.at(indexVec1).R == vector2.at(indexVec2).L){
		                outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			            indexVec1++;
			            indexVec2++;
			            if(indexVec1 == vector1.size()){
			            	STATE = 4;	//vec 1 is empty
			            }
			            else if(indexVec2 == vector2.size()){
			            	STATE = 5;	//vec 2 is empty
			            }
						else{
		                	STATE = 0; //1 = L, 2 = L
						}
		            }
		            else if(vector1.at(indexVec1).R < vector2.at(indexVec2).R){
			            indexVec1++;
		                STATE = 3; //1=L, 2=R
		            }
		            else{//vector1.at(indexVec1).R > vector2.at(indexVec2).R
			            indexVec2++;
		                STATE = 2; //1=R, 2=L
		            }
		            break;//case 1
		        }
		            
		        case 2:{	//1 = R 2 = L
					if(indexVec2 == vector2.size()){
			            outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			            indexVec1++;
						STATE = 5;
					}
					else{
			            if(vector1.at(indexVec1).R == vector2.at(indexVec2).L){
			                indexVec1++;
			            	STATE = 3;
			            }
			            else if(vector1.at(indexVec1).R < vector2.at(indexVec2).L){
			                outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			                indexVec1++;
			                STATE = 0;
			            }
			            else{
			                STATE = 1;
			        	}
					}
		            break;//case 2
		        }
		            
				case 3:{ //1 = L,  2 = R
					if(indexVec1 == vector1.size()){
			            outputVector.at(indexVecOut).R = vector2.at(indexVec2).R;
			            indexVec2++;
						STATE = 4;
					}
					else{
			            if(vector1.at(indexVec1).L == vector2.at(indexVec2).R){
			                indexVec2++;
			            	STATE = 2;
			            }
			            else if(vector1.at(indexVec1).L < vector2.at(indexVec2).R){
			                STATE = 1;
			            }
			            else{
			                outputVector.at(indexVecOut).R = vector2.at(indexVec2).R;
			                indexVec2++;
			                STATE = 0;
			        	}
					}
		            break;//case 3
		        }
		        case 4:{ //vector1 is empty, copy vector2
		            while(indexVec2 < vector2.size()){
						Segment tempDoublet;
						tempDoublet.L = vector2.at(indexVec2).L;
						tempDoublet.R = vector2.at(indexVec2).R;
						outputVector.push_back(tempDoublet);
		                indexVec2++;
		            }
		            STATE = 6; // finished
		            break;//case 4
		        }
		        case 5:{ //vector2 is empty, copy vector1
		            while(indexVec1 < vector1.size()){
						Segment tempDoublet;
						tempDoublet.L = vector1.at(indexVec1).L;
						tempDoublet.R = vector1.at(indexVec1).R;
						outputVector.push_back(tempDoublet);
		                indexVec1++;
		            }
		            STATE = 6; // finished
		            break;//case 5
		        }
		        case 6 :{ //both vectors are empty.
		            finished = true;
		            break;//case 6
				}
			} // end switch
		} // end while
	}// end else
	
	return outputVector;
}

