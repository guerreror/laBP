/*
 *  typedefs.h
 *  InvCoal
 *  Version:  labp_v17
 *
 *
 */
 
 
#ifndef TYPEDEFS_
#define TYPEDEFS_

#include <map>
using std::map;

// Enumerators (used e.g. by the structure Context):

enum Sex {F, M};				// Used for Sex-of-carrier (= SOC) and Sex-of-origin (= SOO)
enum Karyotype {SS, SI, II};	// SS = standard homokaryotype, SI = heterokaryotype, II = inverted homokaryotype
enum Base {A, C, G, T};			// The four DNA bases


// Structures:

struct Segment
//
//	Represents a single segment of chromosome by its left and right map boundaries
//
{
	double L;	// Left boundary
	double R;	// Right boundary
};


struct Context 

{
public:
	unsigned int	pop;	// Population number
	unsigned short inversion; //state of inversion (0=standard)
	
	Context(){}
	Context(const int B,const int C) {pop=B; inversion=C;}
	Context(const Context& other) { pop=other.pop; inversion= other.inversion;}
	
	bool operator<(const Context& other) const    {
		if(pop<other.pop) return true;
		if(pop==other.pop && inversion < other.inversion) return true;
		return false; }
	bool operator==(const Context& other) const    {
		if(pop==other.pop && inversion == other.inversion) return true;
		return false;
	}

};

typedef map<Context, int> cluster_t;


#endif // TYPEDEFS_


// ==================================================================



//Debugging Code
#ifdef DEBUGGER
#define DBG(x) cout << x << endl;
#else
#define DBG(x)
#endif

//Current Debugging so everything isn't so messy Code
#ifdef CURRENTDEBUGGER
#define CDBG(x) cout << x << endl;
#else
#define CDBG(x)
#endif

//Recombination Debug
#ifdef RCDBGR
#define RCDBG(x) cout << x << endl;
#else
#define RCDBG(x)
#endif

//Splitting Segs Debugging Code
#ifdef SRCDBGR
#define SRCDBG(x) cout << x << endl;
#else
#define SRCDBG(x)
#endif



//Query Debug
#ifdef QDBGR
#define QDBG(x) cout << x << endl;
#else
#define QDBG(x)
#endif

//Coelescence Debug
#ifdef COALDBGR
#define COALDBG(x) cout << x << endl;
#else
#define COALDBG(x)
#endif
