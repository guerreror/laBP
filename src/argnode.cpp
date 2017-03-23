/*
 *  argnode.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *  Started by Mark Kirkpatrick on 7/9/08
 *  Edited by Rafael Guerrero 2008-2016
 *
 */
#include <memory>
using std::shared_ptr;


#include <iostream>
	using std::cout;
	using std::endl;
	
#include "chromosome.h"
#include "argnode.h"
#include "typedefs.h"




ARGNode::ARGNode()	// Null constructor.  Used to make an (undefined) ancestor for an ARGNode
{
}

ARGNode::ARGNode(unsigned long dummy)	// dummy constructor.  Used to make a descendant for a terminal ARGNode
{
	nodeNumber=-1;
	time=-1;
	
}

ARGNode::ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr)
//
// Constructor for a terminal node.  The segments carried by this node are
//	stored in descendant.segmentVec[0]
{	
	shared_ptr < ARGNode > first;
	first.reset(new ARGNode(0));
	nodeNumber = nNode;
	time = 0;
	context = chr->getContext();
	Branch branch(first, chr->get_All_SegVec());
	descendant.push_back( branch );
}



ARGNode::ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr, const double t)
//
// Constructor for a recombination node.  Also sets the current node as
//	the ancestor of the chromosome's descendant node.
{
	nodeNumber = nNode;
	time = t;
	context = chr->getContext();
	Branch branch (chr->getDescendant(), chr->get_All_SegVec() );
	descendant.push_back( branch );	
	chr->getDescendant()->addAncestor(this);
}



ARGNode::ARGNode(unsigned long nNode, shared_ptr<Chromosome> chr1, shared_ptr<Chromosome> chr2, const double t)
//
// Constructor for a coalescent node.  Sets the current node as
//	the ancestor of each chromosome's descendant node.
//
{
	nodeNumber = nNode;
	time = t;
	context = chr1->getContext();
	Context context2 = chr2->getContext();
	if(!(context == context2)) std::cerr <<"\nError ARGNode::ARGNode ===> coalescence between carriers of different contexts\n";
	Branch branch1 (chr1->getDescendant(), chr1->get_All_SegVec());
	descendant.push_back( branch1 );	
	Branch branch2 (chr2->getDescendant(), chr2->get_All_SegVec());
	descendant.push_back( branch2 );
	chr1->getDescendant()->addAncestor(this);
	chr2->getDescendant()->addAncestor(this);
}



ARGNode::~ARGNode()	// Destructor
{
	// For debugging:
	//	cout << "  ARGNode::~ARGNode: Deleting node " << nodeNumber << endl;
}



void ARGNode::addAncestor(ARGNode * ancNode)	// Add ancNode as an ancestral node
{
	// For debugging:
	//	cout << "  ARGNode::addAncestor: Adding node " << ancestor[nAncestors]->getNodeNumber() << 
	//			" as ancestor " << nAncestors << " of node " << nodeNumber << endl;

	ancestor.push_back( ancNode );
	if(ancNode->getTime() < time)
        std::cerr << endl << ">>>WARNING: ARGNode::addAncestor: ancestor at time " << ancNode->getTime() <<
						" is younger than its descendant node at time " << time << endl << endl;
}



void ARGNode::addDescendant(shared_ptr < ARGNode >descNode, vector<Segment> segVec)
{
//
// Adds descNode as a descendant ARG node, and segVec as its corresponding segment array.
//	Also sets the current node as the ancestor of that descendant.
//

	// For debugging:
	//	cout << "  ARGNode::addDescendant: Adding node "<< descendant[nAncestors]->getNodeNumber() << 
	//			" as descendant " << nDescendants << " of node " << nodeNumber << endl;

	if(descNode->getTime() > time)
		std::cerr << endl << ">>>WARNING: ARGNode::addDescendant: descendant at time " << descNode->getTime() <<
						" is older than its ancestor node at time " << time << endl << endl;
						
	Branch branch (descNode, segVec);
	descendant.push_back( branch );	
	descNode->addAncestor(this);				// Set the ancestor of descNode equal to this node
}
	


unsigned long ARGNode::getNodeNumber()	// Return the identifier number of this node
{
	return(nodeNumber);
}



double ARGNode::getTime()	// Return the time (= age) of this node
{
	return(time);
}



Context ARGNode::getContext()	// Return the context this node
{
	return(context);
}



unsigned long ARGNode::getNAncestors()	// Return the number of ancestral nodes from this node
{
	return(ancestor.size());
}



ARGNode* ARGNode::getAncestor(unsigned long i)	// Get ancestral node i
{
	// For debugging:
	//	cout << "  ARGNode::getAncestor: Getting ancestor[" << i << "] of node " << nodeNumber << endl;

	if(i>ancestor.size()-1)		// Could handle this error more gracefully?
	{
		cout << endl << ">>>ERROR: ARGNode::getAncestor called for ancestor[" << i << "], but this node has only ";
		cout << ancestor.size() << " ancestors" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(ancestor[i]);
}



unsigned long ARGNode::getNDescendants()	// Return the number of descendant nodes from this node
{
	return(descendant.size());
}



shared_ptr<ARGNode> ARGNode::getDescendantNode(unsigned long i)	// Get descendant ARG node i
{
	// For debugging:
	//	cout << "   ARGNode::getDescendant:  Getting descendant[" << i << "]..." << endl;

	if(i > descendant.size())	// Can we handle this error more gracefully?
	{
		cout << endl << ">>>ERROR: ARGNode::getDescendant called for descendant[" << i << "], but this node has only ";
		cout << descendant.size() << " descendants" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(descendant.at(i).node);
}



vector< Segment > ARGNode::getDescendantSegmentVector(unsigned long i)
//
// Get the vector of chromosome segments corresponding to the branch leading to descendant i
//
{
	// For debugging:
	//	cout << "     ARGNode::getDescendantSegmentList:  Getting descendant segment list " << i << "..." << endl;

	if(i > descendant.size()-1)		// Is there a better way to handle this error?
	{
		cout << endl << ">>>ERROR: ARGNode::getDescendantSegmentList called for descendantSegmentList[" << i << "], but this node has only ";
		cout << descendant.size() << " descendants" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(descendant.at(i).segmentVec);
}



void ARGNode::outputNode()	// Prints out the node's data
{
	cout << " ARGNode " << nodeNumber << ":" << endl;
	cout << "  time = " << time << endl;
	cout << "  nAncestors = " << ancestor.size() << endl;
	
	if(ancestor.size() > 0)
	{
		cout << "   Ancestal node number(s): ";
		cout << ancestor[0]->getNodeNumber();
		
		for(unsigned long i=1; i<ancestor.size(); i++)
			cout << ", " << ancestor[i]->getNodeNumber();
		cout << endl;
	}
	cout << "  nDescendants = " << descendant.size() << endl;
	if(descendant.size() > 0)
		for(unsigned long i=0; i<descendant.size(); i++)
		{
			
			cout << "   descendant[" << i << "] = node " << descendant.at(i).node->getNodeNumber()  <<  endl;
			cout << "     Segment(s) descending to that node:";
			for(unsigned long j=0; j<descendant.at(i).segmentVec.size(); j++)
				cout << "  (" << descendant.at(i).segmentVec.at(j).L << ", " << descendant.at(i).segmentVec.at(j).R << ")";
			cout << endl;
			}
		
}



void ARGNode::outputARG()	// Prints out the ARG descending from this node
{
	cout << " ARGNode " << nodeNumber << ":" << endl;
	cout << "  time = " << time << endl;
	cout << "  nAncestors = " << ancestor.size() << endl;
	
	if(ancestor.size() > 0)
	{
		cout << "   Ancestal node number(s): ";
		cout << ancestor[0]->getNodeNumber();
		
		for(unsigned long i=1; i<ancestor.size(); i++)
			cout << ", " << ancestor[i]->getNodeNumber();
		cout << endl;
	}
	cout << "  nDescendants = " << descendant.size() << endl;
	if(descendant.size() > 0)
		for(unsigned long i=0; i<descendant.size(); i++)
		{
			
			cout << "   descendant[" << i << "] = node " << descendant.at(i).node->getNodeNumber()  <<  endl;
			cout << "     Segment(s) descending to that node:";
			for(unsigned long j=0; j<descendant.at(i).segmentVec.size(); j++)
				cout << "  (" << descendant.at(i).segmentVec.at(j).L << ", " << descendant.at(i).segmentVec.at(j).R << ")";
			cout << endl;
			
		}
	cout << endl;
	for(unsigned long i=0; i<descendant.size(); i++)
	{
		descendant.at(i).node->outputARG();
	}
}

Branch::Branch(shared_ptr<ARGNode> n, vector<Segment> segVec){
	node= n;
	segmentVec= segVec;
}

Branch::~Branch(){
	
}

Branch::Branch(const Branch& br){
	node=br.node;
	segmentVec= br.segmentVec;
	
}	

void Branch::operator =(const Branch& br){
	node=br.node;
	segmentVec= br.segmentVec;
	
}


