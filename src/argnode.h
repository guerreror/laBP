/*
 *  argnode.h
 *  InvCoal
 *  Version:  labp_v17
 *
 *  Ancestral Recombination Graph class started by Mark Kirkpatrick on 7/9/08.
 *
 *  Adapted for inverted chromosomes by Rafael Guerrero 2008-2016

 *	   An ARGNode is the basic unit of the ancestral recombination graph (ARG).  The call:
 *
 *			ARGNode node0 = ARGNode(nodeID, ctxt, segVec);
 *
 *	makes a new terminal node node0 (at time = 0) that has identifier nodeID (an int), context ctxt,
 *  and vector of segments segVec.  (The nodeID is useful for output, debugging, etc., but has no
 *	internal significance;  no error occurs if multiple nodes have the same nodeID.)  The
 *  call:
 *
 *			ARGNode node1 = ARGNode(nodeID, ctxt, t);
 *
 *	makes node1 as a new root or internal ARGNode node that has descendant ARG node descNode.
 * 
 *	   When an ARGNode is made, no ancestor or descendant nodes are set.  To add the ARGNode node0
 *  as a descendant of node1 (for example), use:
 *
 *			node1.addDescendant(&node0, segVec);
 *
 *	where segVec is the vector of segments carried by the chromosome that links those two nodes.  That  
 *	function call also automatically sets node1 as an ancestor node of node0.  Alternatively, the
 *	ancestor of node0 can be set explicitly using
 *
 *			node0.addAncestor(&node1);
 *
 *		In the case of a terminal node, descendant[0].segmentVec has the segments carried.
 *
*/


#ifndef ARGNODE_
#define ARGNODE_

#include <vector>
	using std::vector;
	
#include <memory>
using std::shared_ptr;
	
#include "typedefs.h"
#include "chromosome.h"

// Forward declarations:
	class Chromosome;
	class ARGNode;


struct Branch
//
// Structure used to represent descendant branches in ARGNode
//
{
	shared_ptr < ARGNode > node;	// Pointer to the descendant ARG node
	vector< Segment > segmentVec;	// Vector of segments carried by the chromosome along this branch
	Branch(shared_ptr < ARGNode > n, vector< Segment > segVec);
	~Branch();
	Branch( const Branch& br );			// Copy constructor (copies the branch)
	void operator = ( const Branch& br );	// Eq
};



class ARGNode
{
  private:
	unsigned long nodeNumber;		// Identifier number of this node
	double time;					// Time at which this node occurred
	Context context;				// Context at this node
	vector<ARGNode*> ancestor;		// Vector of pointers to ancestral ARG nodes
	vector<Branch> descendant;		// Vector of branches descending from this ARG node;  these structures
									//	have pointers to the descendant ARG nodes and the chromosome segments
									//	carried along those branches of the ARG
													
  public:
// Constructors and destructor:
	ARGNode();										// Null constructor;  used to make (empty) ancestor node
	ARGNode(unsigned long i);
	ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr);		// Constructor for a terminal node
	ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr,
			double t);								// Constructor for a recombination node
	ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr1,
			shared_ptr<Chromosome> chr2, double t);	// Constructor for a coalescent node
	~ARGNode();										// Destructor
// Adding ancestors and descendants to a node:
	void addAncestor(ARGNode * ancNode);			// Adds ancNode as an ancestral ARG node
	void addDescendant(shared_ptr < ARGNode > descNode,			// Adds descNode as a descendant node and adds segList as its corresponding 
						vector<Segment> segList);	//	vector of segments;  sets current node as ancestor of the new descendant
// Functions that get values:
	unsigned long getNodeNumber();							// Get the identifier number of this node
	double getTime();									// Get the time (= age) of this node
	Context getContext();							// Get the context for this node
	unsigned long getNAncestors();							// Get the number of ancestral nodes of this node
	ARGNode * getAncestor(unsigned long i);					// Get pointer to ancestral node i
	unsigned long getNDescendants();							// Get the number of descendant nodes from this node
	shared_ptr < ARGNode > getDescendantNode(unsigned long i);				// Get pointer to descendant node i
	vector<Segment> getDescendantSegmentVector(unsigned long i);	// Get descendant segment vector i
// Output data for this node:
	void outputNode();								// Prints out the node's data
	void outputARG();								// Prints out the ARG descending from this node 
};


#endif //ARGNODE_
