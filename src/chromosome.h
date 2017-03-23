/*
 * Chromosome.h
 *  InvCoal
 *  Version:  labp_v17
 *
 *
 */


#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

// Includes for STL:
	#include <vector>
		using std::vector;

    #include <memory>
		using std::shared_ptr;
		using std::unique_ptr;

// Includes for our files:
	#include "typedefs.h"
	#include "argnode.h"


// Forward declaration:
	class ARGNode;


class Chromosome 
{
	
  private:
	struct ChromosomeData;								// Contains all of the chromosome's private data
	unique_ptr<ChromosomeData> chromData;				// Pointer to the private data
	vector<Segment> combineDoubletVectors(vector<Segment> vector1, // Combines two segment vectors into one
							vector<Segment> vector2);

  public:
	// Constructors and destructor:
	Chromosome();										// Empty constructor				
	Chromosome( Context ctxt,							// The typical constructor
				vector<Segment> chromSegments,			
				shared_ptr < ARGNode > desc);
	Chromosome( Context ctxt,							//constructor for chromosome with inversion
				vector<Segment> chromSegments,
				Segment invRange,
			   shared_ptr < ARGNode > desc );
	~Chromosome();										// Destructor
	
	// Copy functions:
	Chromosome( Chromosome& chr );			// Copy constructor (copies the chromosome)
	void operator = ( Chromosome& chr );	// Equals operator;  has same effect as the copy constructor	
	shared_ptr<Chromosome> copy_values(Segment& inv_range);
	// Set and get data for this chromosome:
	void setPopulation(unsigned short pop);           // Sets population number
	void setInvSegs(vector<Segment> segs);
	void setStdSegs(vector<Segment> segs);
	void setInversion(unsigned short inv);
	void setDescendant(shared_ptr < ARGNode > desc);		// Set the pointer to the descendant ARG node
	void setContext(Context c);
	void clearInvSegs();

	Context getContext();					// Returns whole context
	unsigned short getPopulation();					// Returns population number
	unsigned short getInv(); 
	shared_ptr < ARGNode > getDescendant();	// Returns pointer to the descendant ARG node
	vector<Segment> get_All_SegVec();		// Returns vector of chromosome segments
	vector<Segment> get_Inv_Segs();
	vector<Segment> get_Std_Segs();
	double getHomoLength(Segment& invRange);					// Returns length of the chromosome	available for homokaryotypic recombination
	double getHeteroLength(Segment& invRange);
	Segment getHomoRange(Segment& invRange); 
	// Recombination and coalescent functions:
	vector<Segment> translateSegs(const Segment& invRange);
	double calcBreakpoint(Segment& invRange,bool hetero);
	bool isEmpty();
    bool carriesSite(double x);
	void merge(shared_ptr <Chromosome> chrom, shared_ptr < ARGNode > descNode);				// Adds segments of chrom to this chromosome
	shared_ptr<Chromosome> recombine(double breakpoint, 
									 shared_ptr< ARGNode > descNode, 
									 Context newContext, 
									 Segment& invRange
									 );
	
	shared_ptr<Chromosome> doubrecombine(double breakpoint1, 
										 double breakpoint2, 
										 shared_ptr< ARGNode > descNode, 
										 Context newContext, 
										 Segment& invRange
										 );
	
};	

vector<Segment> splitSegments(vector<Segment>& old, double breakpoint, bool inverted);

struct Chromosome::ChromosomeData{
	//
	// This structure contains all the data for a chromosome
	//	
	Context context;						// Contains the chromomsome's SOO, SOC, population, etc.
	vector<Segment> chromosomeSegments;		// Vector of chromosome segments, each of which is a pair of floats
	vector<Segment> invertedSegments;	
	shared_ptr< ARGNode > descendant;					// Pointer to the descendant ARG node of this chromosome
};

double tr_pt(double p, const Segment& invRange);

#endif /*CHROMOSOME_H_*/
