/*
 *  sitenode.h
 *
 *  Created by Mark Kirkpatrick on 7/28/08, modified 4-IX-08 so arrays are allocated dynamically using vector<>.
 *
 *	The SiteNode class contains the data needed to represent a node in a gene tree pertaining
 *	to a site on a chromosome whose evolution is represented by an ancestral recombination graph (ARG).
 *	Definitions for the members and methods are given below.
 *
 *  The call:
 *				geneTree = SiteNode(sitePos, &argNd, NULL)
 *
 *	returns the SiteNode for the gene tree for the site sitePos using as root the ARGNode argNd.  The 
 *	call:
 *
 *				geneTree = SiteNode(sitePos, &argNd, &siteNd)
 *
 *	returns the SiteNode that is the root for the subtree descending from argNd and setting its ancestral
 *	site node to the SiteNode siteNd.
 *
 *  baseState is an array of mutations that lie between this node and the SiteNodes ancestral to it.
 *	it.  mutationList[0] is the current base at this site, even if there were no mutations along
 *	the ancestral branch.  mutationList[1], mutationList[2], ... list mutations (bases) that occurred 
 *	along that branch in reverse chronological order.
 *
 */

#ifndef SITENODE_
#define SITENODE_

#include "typedefs.h"
#include "argnode.h"

struct snpHit;

class SiteNode
{
private:
	double sitePosition;				// The position of the site, typically given in morgans
	unsigned long nodeNumber;					// The identifier number of this node;  agrees with the nodeNumber for the corresponding ARGNode
	double time;						// Time (age) of the node
	double branchL;
    double branchL_Informative;
    Context context;				// Context of the chromosome at this node
	Context original_context;
	vector<shared_ptr < SiteNode > > descendant;	// Vector of pointers to descendant nodes;  descendant.size() returns the number of descendants.
	vector<Base> baseState;			// Array with the states on the branch ancestral to this SiteNode.  baseState[0] = baseState.front() is
									//	the base state of this node's ancestral node, baseState.back() is the base state at this node,
									//	baseState.size() is the number of base states in this vector, baseState.size() - 1 is the 
									//	number of mutations that happened between this node and its ancestor.
	//
    // Functions that work with mutation:
	vector<Base> mutate_jc(Base base0, double mu, 
							double time);			// base0 is the base state at the ancestral node, mu the mutation rate,
													//	and time the number of generations between this node and its ancestral node
    void putBaseStates(vector<Base> baseVec);		// Sets baseState (the vector of ancestral base states) equal to baseVec

public:
	// Constructors and destructor:
	SiteNode();											// Null constructor
	SiteNode(double sitePos, shared_ptr < ARGNode > argRoot);			// Typical constructor
	~SiteNode();										// Destructor
	//
	// Function that puts mutations along branches:
	void mutateDescendants(Base base0, double mu);		// Put mutations on all branches of the gene tree descending from this node,
														//	given this node has base base0, a mutation rate mu

  	// Functions that return data about this node:
	unsigned long getNodeNumber();								// Returns the node number
	double getTime();									// Returns the time (age) of the node
	double getBranchLength();
    double getBranchLength_Informative();
    Context getContext();								// Returns the context of the node
	unsigned int getNDescendants();								// Return the number of descendants
	shared_ptr < SiteNode> getDescendant(unsigned int i);					// Return a pointer to descendant i
	Base getBaseState();								// Get the base state at this node;  equivalent to getBaseState.back()
	Base getBaseState(unsigned int i);							// Get base state i along ancestral branch

	//
	// Output data for this node and/or the whole gene tree:
	void outputNode();									// Output data for this node only
	void outputTree();									// Output the gene tree descending from this node
    void outputShortTree();
    vector<unsigned>   outputSNPs(vector<unsigned>  s);
    double getTotalLength(double runtot);
    double getTotalLength_Informative(double runtot);
    double getInvLength(double runtot);
    void setBranch(double ancTime);
    void setBranch_informative(double ancTime);
    int calcBranchLengths(int nCall);
    int calcBranchLengths_informative(int nCall);
    snpHit getSNPhit(double target, snpHit x);

};

struct snpHit{
    double runningTot;
    bool found;
    shared_ptr<SiteNode> chosen;
};
#endif //SITENODE_
