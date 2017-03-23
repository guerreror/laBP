/*
 *  siteNode.cpp
 *
 *  Created by Mark Kirkpatrick on 7/28/08.
 *
 *	See "siteNode.h" for documentation.
 *
 */

#include <iostream>
	using std::cout;
	using std::endl;
	
#include "sitenode.h"
#include "argnode.h"
#include "ran_mk.h"


static bool siteQuery(double sitePos, vector<Segment> segList)
//
//	Returns TRUE if the chromosome at the ARG node *argNode carries site sitePos, otherwise returns FALSE
//
{
    for(int i=0; i<segList.size(); ++i)
    {
        // For debugging:
        QDBG( "      siteNode::siteQuery:  "<<i<<" " << segList[i].L << " < " <<
             sitePos << " < " << segList[i].R << " ? " )
        
        if( (segList.at(i).L < sitePos) && (sitePos < segList.at(i).R) )
        { QDBG("					YES")
            return(true);
        }
        else if( segList.at(i).L == sitePos )
        { QDBG("					YES")
            return(true);
        }
        
    }
    
    return(false);
} // end siteQuery

static shared_ptr < ARGNode > getNextSiteNode(double sitePos, shared_ptr < ARGNode > argRoot)
//
//	Descends the gene tree for the site sitePos starting at the shared_ptr < ARGNode > argRoot until the next split
//	or a terminal node is reached, then returns that ARGNode.
//
{
    // Count the number and make a list of the descendant ARGNodes from *argRoot that carry our site:
    //
    int uniqueDesc=0, nSiteDesc = 0, nARGDesc = argRoot->getNDescendants();
    for(int i = 0; i < nARGDesc; i++)
        if(siteQuery(sitePos, argRoot->getDescendantSegmentVector(i)) )
        {
            nSiteDesc++;
            uniqueDesc = i;	// If there is only one descendant of this ARGNode that carries this site, uniqueDesc
            //  will be the index number of that descendant in the descendant array of argRoot.
        }
    QDBG("Query says "<<nSiteDesc<<" and unique is "<<uniqueDesc<<'\n')
    // Process the ARG node depending on the number of descendants:
    //
    if(nSiteDesc > 1)	{// *argRoot is a split
        QDBG("nSite is >1, returning argRoot\n")
        return(argRoot);
    }
    
    else if(nSiteDesc == 1 && argRoot->getTime() != 0)	// *argRoot is not at the root for the gene tree for this site, so we proceed
    {					//	to the descendant ARGNode that carries this site
        QDBG("nSite is ==1, recursive step...\n")
        return(getNextSiteNode(sitePos, argRoot->getDescendantNode(uniqueDesc)));
    }
    // If we reach here, nSiteDesc == 0, and this is a terminal node:
    //
    QDBG("nSite is ==0, is this a terminal node?\n")
    return(argRoot);
    
} // end getNextSiteNode


SiteNode::SiteNode()	// Null constructor
{
}



SiteNode::SiteNode(double sitePos, shared_ptr < ARGNode > argRoot)		// The typical constructor
//
// Constructs the gene tree for a site at chromosome position sitePos that descends from the ARG node *argRoot
// and that has ancestral site node *ancSiteNode.  To construct the gene tree from the root of its ARG,
// call:
//			SiteNode(sitePos, *argRoot)
//
//
//
{
	// For debugging:
	QDBG("\n   Making a SiteNode at node " << argRoot->getNodeNumber() )

	// Descend the gene tree until we reach the next split at our site:
	// 
	shared_ptr < ARGNode > argSiteRoot = getNextSiteNode(sitePos, argRoot);	// Find the ARGNode that corresponds to the root for the gene tree at this site.
	
	// For debugging:
	//	if(argSiteRoot != argRoot)
	//		cout << "    (Skipped to node " << argSiteRoot->getNodeNumber() << ")" << endl;

	// Set the basic private data:
	//
	sitePosition = sitePos;
	nodeNumber = argSiteRoot->getNodeNumber();
	time = argSiteRoot->getTime();
	context = argSiteRoot->getContext();

	// Set the descendant nodes:
	//
	unsigned long argSiteRootNDesc = argSiteRoot->getNDescendants();
	unsigned long i;
	for(i = 0; i < argSiteRootNDesc; i++)
		if(siteQuery(sitePos, argSiteRoot->getDescendantSegmentVector(i)) )				// descendant[i] of *argSiteRoot is a carrier for sitePosition, 
		{
			shared_ptr<SiteNode> d;
			d.reset(new SiteNode(sitePos, argSiteRoot->getDescendantNode(i)));
			descendant.push_back(d);															//	so set that node as a descendant of this node and make a 
		}																						//	new SiteNode for it, which then recursively descends down that 
																								//	branch of the gene tree.
} // end SiteNode::SiteNode



SiteNode::~SiteNode()		// Destructor
{
	// For debugging:
	//	cout << "  SiteNode::~SiteNode called for node " << nodeNumber << endl;
}



vector<Base> SiteNode::mutate_jc(Base base0, double mu, double t)
//
//	Returns a vector with the base states along a descendant branch from a node
//	with base state base0, given a mutation rate mu and the number of generations
//	between the nodes is t.
//
//	Uses the Jukes-Cantor model:  changes between all pairs of bases are equally likely.
//
{
	vector<Base> baseVec(1, base0);		// Make a vector of bases and set the first element to base0 
	double timeLeft = t;
	
	do{
		double waitTime = randexp(mu);	// Waiting time until the next mutation, assumed independent of the current base state
		
		if(timeLeft < waitTime)			// No more mutation, return
			return baseVec;
			
		else {							// Mutation happens
			int newBase = ( baseVec.back() + randint(1, 3) ) % 4;	// All changes assumed equally likely
			baseVec.push_back(Base(newBase));						// Converts newBase, which is an int, to a Base via a type case
		} 
			
		timeLeft -= waitTime;
		
	} while(timeLeft > 0);	// Normally a value is returned before we get here
	
	return baseVec;
}



void SiteNode::mutateDescendants(Base base0, double mu)
//
//	Puts mutations on the descendant branches of this node (that is, sets the values of the vector baseState
//	  in the descendant nodes), given base0 as the base state at the ancestral node, mu as the mutation rate
//
{
	static unsigned long nCalls;	// Equals 0 i.f.f. this is the first time the rountine is called
	if(nCalls == 0) {	// ... then this node is the root of the gene tree, so set its base state to base0
		baseState.push_back(base0);
		nCalls++;
	}	
	
	for(unsigned long iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
	{
		vector<Base> baseVec = mutate_jc(base0, mu, time - descendant[iDesc]->getTime());	// Make vector of base states along the branch to descendant iDesc
		descendant[iDesc]->putBaseStates(baseVec);											// Set the vector baseState of descendant[iDesc] equal to that vector
		descendant[iDesc]->mutateDescendants(descendant[iDesc]->getBaseState(), mu);		// Descend gene tree recursively to that descendant and repeat
	}
}



void SiteNode::putBaseStates(vector<Base> baseVec)
//
// Sets baseState (the vector of ancestral base states) equal to baseVec
//
{
	baseState = baseVec;
}



unsigned long SiteNode::getNodeNumber()
{
	return(nodeNumber);
}



double SiteNode::getTime()
{
	return(time);
}


double SiteNode::getBranchLength()
{
    return(branchL);
}


Context SiteNode::getContext()
{
	return(context);
}



unsigned int SiteNode::getNDescendants()
{
	return(descendant.size());
}



shared_ptr< SiteNode> SiteNode::getDescendant(unsigned int i)
{
	return(descendant[i]);
}



Base SiteNode::getBaseState()
{
	return(baseState.back());
}



Base SiteNode::getBaseState(unsigned int i)
{
	return(baseState.at(i));
}


void SiteNode::outputNode()
{
	cout << " SiteNode " << nodeNumber << " for site " << sitePosition << endl;
	cout << "   time = " << time << endl;
	
	if(descendant.size() == 0)
		cout << "   Terminal node." << endl;
	else
	{
		cout << "   nDescendants = " << descendant.size() << endl;
		cout << "     Node(s):  ";
		for(unsigned long i = 0; i < descendant.size(); i++)
			cout << descendant[i]->getNodeNumber() << "  ";
		cout << endl;
	}
	cout << "   number of base states on ancestral branch = " << baseState.size() << endl << "   ";
	
	for(unsigned long i=0; i < baseState.size(); i++)
		cout << "  " << baseState[i];
	if(baseState.size() > 0)
		cout << endl;
	cout << endl;
}


void SiteNode::outputTree()
{
	cout << " SiteNode " << nodeNumber << " for site " << sitePosition << endl;
	cout << "   time = " << time << endl;
	
	if(descendant.size() == 0)
		cout << "   Terminal node" << endl;
	else
	{
		cout << "   nDescendants = " << descendant.size() << endl;
		cout << "     Node(s):  ";
		for(unsigned long i = 0; i < descendant.size(); i++)
			cout << descendant[i]->getNodeNumber() << "  ";
		cout << endl;
	}
	cout << "   number of base states on ancestral branch = " << baseState.size() << endl << "   ";
	
	for(unsigned long i=0; i < baseState.size(); i++)
		cout << "  " << baseState[i];
	if(baseState.size() > 0)
		cout << endl;
	cout << endl;
	for(unsigned long i=0; i < descendant.size(); i++)
		descendant[i]->outputTree();
}

void SiteNode::outputShortTree()
{
    static int tab=0;
    if(time!=-1){
        cout << "S_" << nodeNumber<<"(" <<time<<")[" << context.inversion <<"] ";
        
        if(time != 0)
        {
            //cout << "   nDescendants = " << descendant.size() << endl;
            cout << "TO: ";
            for(unsigned int i = 0; i < descendant.size(); i++)
                cout << descendant[i]->getNodeNumber() << " ("<<descendant[i]->getBranchLength()<<") ";
            cout << endl;
            
            for(unsigned int i=0; i < descendant.size(); i++){
                tab++;
                for (int x=0; x<tab; ++x) cout<<" ";
                descendant[i]->outputShortTree();
                tab--;
            }
            
        }
        else{cout<<endl;}
        
    }
}

void SiteNode::setBranch(double ancTime){
    if(time!=-1) branchL= ancTime - time;
    else branchL=0;
    
    if(branchL<0)cout<< "Error in setBranch:: setting a negative length\n";
}

void SiteNode::setBranch_informative(double ancTime){
    if(time > 1) branchL_Informative= ancTime - time;
    else branchL_Informative=0;
    
    if(branchL_Informative<0)cout<< "Error in setBranch:: setting a negative length\n";
}

int SiteNode::calcBranchLengths_informative(int nCall){
    
    int nCalls = nCall;
    
    if(nCalls == 0) {
        branchL=0;
        nCalls++;
    }
    
    for(unsigned int iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
    {
        descendant[iDesc]->setBranch_informative(time); //NOTICE THAT IT USES THE INFORMATIVE, BIASED LENGTH!!
        descendant[iDesc]->calcBranchLengths_informative(nCalls);
    }
    return nCalls;
}

int SiteNode::calcBranchLengths(int nCall){
    
    int nCalls = nCall;
    
    if(nCalls == 0) {
        branchL=0;
        nCalls++;
    }
    
    for(unsigned int iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
    {
        descendant[iDesc]->setBranch(time);
        descendant[iDesc]->calcBranchLengths(nCalls);
    }
    return nCalls;
}

double SiteNode::getTotalLength(double runtot){
    
    double runTotal =runtot;
    // cout<<"["<<nodeNumber<<"] "<<runTotal<<'\n';
    for(unsigned int iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
    {
        runTotal+=	descendant[iDesc]->branchL;
        runTotal = descendant[iDesc]->getTotalLength(runTotal);
    }
    
    return runTotal;
}

double SiteNode::getTotalLength_Informative(double runtot){
    
    double runTotal =runtot;
    // cout<<"["<<nodeNumber<<"] "<<runTotal<<'\n';
    for(unsigned int iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
    {
        runTotal+=	descendant[iDesc]->branchL_Informative;
        runTotal = descendant[iDesc]->getTotalLength_Informative(runTotal);
    }
    
    return runTotal;
}

double SiteNode::getInvLength(double runtot){
    
    double runTotal =runtot;
    // cout<<"["<<nodeNumber<<"] "<<runTotal<<'\n';
    for(unsigned int iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
    {
        if(context.inversion==1 ) runTotal+=	descendant[iDesc]->branchL;
        runTotal = descendant[iDesc]->getInvLength(runTotal);
    }
    
    return runTotal;
}


snpHit SiteNode::getSNPhit(double target, snpHit x){
    snpHit c = x;
    
    for(unsigned int i=0; i < descendant.size(); i++){
        if(!c.found){
            c.runningTot+=descendant[i]->getBranchLength();
            
            if(c.runningTot > target) {
                c.found=true;
                c.chosen = descendant[i];
                //cout<<c.runningTot<<" is good, returning node "<<descendant[i]->getNodeNumber()<<"\n";
                return c;
            }
            else{
                // cout<<c.runningTot<<", not yet.\n";
                c = descendant[i]->getSNPhit(target, c);
            }
        }
    }
    return c;
}


vector<unsigned>  SiteNode::outputSNPs(vector<unsigned>  s)
{
    // this function returns 1 for all terminal descendants of a SiteNode
    // It's supposed to be used only on the SiteNode output from getSNPhit in order to get polymorphic loci
    
    vector<unsigned> sample = s;
    if(time == 0)
    {
        sample[nodeNumber]=1;
    }
    for(unsigned int i=0; i < descendant.size(); i++){
        sample = descendant[i]->outputSNPs(sample);
    }
    
    return sample;
}

