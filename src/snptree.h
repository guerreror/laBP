/**
 *  snptree.h
 *  InvCoal
 *  Version:  labp_v17
 *
 *  Created by Rafael Guerrero on 11/8/13.
 *  The SNPtree class holds a gene tree (SiteNode class) and some descriptors (length, heterozygosity, allele states). 
 *  I created it to complement SiteNode, be able to play around with gene trees and make it easier to calculate statistics.
*/

#ifndef __snptree__
#define __snptree__

#include <iostream>
#include <vector>
using std::vector;
#include <algorithm> 
using std::sort;
using std::reverse;
#include <numeric>
using std::accumulate;


#include "sitenode.h"

struct SNPtree {
    SiteNode tree; // gene tree
    double totalLength; //sum of all branches in tree
    double length_Informative; //tree length ignoring terminal branches, useful when interested in LD (singleton SNPs are not informative)
    int nCarriers;
    vector<unsigned> SNPvalues; //allele states for carriers (in present).
    double pi; //heterozygosite
    
    SNPtree(SiteNode geneTree, int nCarr);
    SNPtree(SiteNode geneTree, int nCarr, double theta);
    ~SNPtree(){};
    
};

//Constructor for segregating site (i.e., a SNP will always be dropped somewhere in the tree
SNPtree::SNPtree(SiteNode geneTree, int nCarr){
    tree = geneTree;
    tree.calcBranchLengths(0);
    tree.calcBranchLengths_informative(0);
   // tree.outputShortTree();
    
    nCarriers=nCarr;
    if(nCarriers >2) {length_Informative= tree.getTotalLength_Informative(0);}
    else {
        //std::cerr<< "Warning: Informative sites length doesn't work for 2 carriers\n";
        length_Informative=0;
    }
    totalLength=tree.getTotalLength(0);
    
    vector<unsigned> zeroVec(nCarriers, 0);
   
    double target=randreal(0, totalLength);
    
    shared_ptr< SiteNode> null; null.reset( new SiteNode());
    snpHit dum ={0,false,null};
    shared_ptr< SiteNode> startMutant = tree.getSNPhit(target, dum).chosen;

    SNPvalues = startMutant->outputSNPs(zeroVec);
    
};

//Constructor for stochastic mutation. SNPs generated on tree only at rate theta
SNPtree::SNPtree(SiteNode geneTree, int nCarr, double theta){
    tree = geneTree;
    tree.calcBranchLengths(0);
    tree.calcBranchLengths_informative(0);
    
    nCarriers=nCarr;
    if(nCarriers >2) {length_Informative= tree.getTotalLength_Informative(0);}
    else {
        //std::cerr<< "Warning: Informative sites length doesn't work for 2 carriers\n";
        length_Informative=0;
    }
    totalLength=tree.getTotalLength(0);
    
    if(totalLength==0){std::cerr <<"Error: In SNPtree, total gene tree length is zero.\n";exit(1);}
    
    double target=randreal(0, totalLength);
    double waitingT = 0;
    if(theta > 0) waitingT = randexp(theta);
    
    vector<unsigned> zeroVec(nCarriers, 0);

    if (waitingT < totalLength){
        shared_ptr< SiteNode> null; null.reset( new SiteNode());
        snpHit dum ={0,false,null};
        shared_ptr< SiteNode> startMutant = tree.getSNPhit(target, dum).chosen;
        SNPvalues = startMutant->outputSNPs(zeroVec);
    }
    else{
        SNPvalues = zeroVec;
    }
};

// goes through a collection of SNPtrees (i.e., sites with potentially different gene trees), gets a haplotype per carrier, determines how many copies of that haplotype
// Useful for allele frequency spectrum
vector<int> alleleCount(vector< SNPtree > trees){
    
    vector<vector<int> > haplotypes;
    haplotypes.resize(trees[0].nCarriers);
    
    for(int i=0; i< trees[0].nCarriers;++i){
        for(int k=0; k<trees.size() ;++k){
            haplotypes[i].push_back(trees[k].SNPvalues[i]);
        }
    }
    
    unsigned long int nCarriers= haplotypes.size();
    
    vector<int> tally(nCarriers,0); // We'll fill this with the counts of unique haplotypes
    vector<vector<int> > unique; //The list of unique haplotypes
    unique.push_back(haplotypes[0]);
    tally[0]=1; //Arbitrarily, the first haplotype has at least one copy (the first carrier
    
    for(int x=1; x<nCarriers; ++x){
        bool match =false;
        unsigned long int match_id = unique.size();
        
        for (int y=0; y< unique.size() ;++y){
            if(unique[y] == haplotypes[x]){
                match=true;
                match_id=y;
        //        cout<<"Unique["<<y<<"] matches haplo["<<x<<"]\n";
            }
        }
        
        if(!match){unique.push_back(haplotypes[x]);}
        tally[match_id]++;
    }
    
    sort(tally.begin(), tally.end());
    reverse(tally.begin(), tally.end());
  
    vector<int> noZeroes;
    for(int i =0 ; i< tally.size(); ++i) {if(tally[i]!=0) noZeroes.push_back(tally[i]);}
    
  // for(int i=0; i<noZeroes.size();  ++i) {cout<<noZeroes[i]<<" ";}cout<<'\n';

    return noZeroes;
};

#endif /* defined(__snptree__) */
