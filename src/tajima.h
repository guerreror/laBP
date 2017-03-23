/**
 *  tajima.h
 *  InvCoal
 *  Version:  labp_v17
 * 
 *
 *  Functions for population genetic statistics
 *  The functions for tajima's D were written by Monty Slatkin
 *  Other functions (heterozygosity, Dxy, Fst) were written by Rafael Guerrero 
*/

#ifndef tajima_h
#define tajima_h

#include <stdio.h>
#include <math.h>

/*
 ************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
 and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
 D (a1, a2, b1, b2, c1, c2, e1, e2).
 *************************************************************************************************
*/

double a1f(int nsam)
{
    double a1;
    int i;
    a1 = 0.0;
    for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
    return (a1);
}
double a2f(int nsam)
{
    double a2;
    int i;
    a2 = 0.0;
    for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
    return (a2);
}
double b1f(int nsam)
{
    double b1;
    b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
    return (b1);
}
double b2f(int nsam)
{
    double b2;
    b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
    return (b2);
}
double e1f(double a1, double c1)
{
    double e1;
    e1 = c1/a1;
    return (e1);
}
double e2f(double a1, double a2, double c2)
{
    double e2;
    e2 = c2/((a1*a1)+a2);
    return (e2);
}
double c1f(double a1, double b1)
{
    double c1;
    c1 = b1 - (1/a1);
    return (c1);
}
double c2f(int nsam, double a1, double a2, double b2)
{
    double c2;
    c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
    return (c2);
}
////////////////////////////
double tajd(int nsam, int segsites, double sumk)
{
    
    double  a1, a2, b1, b2, c1, c2, e1, e2;
    
    if( segsites == 0 ) return( 0.0) ;
    
    a1 = a1f(nsam);
    a2 = a2f(nsam);
    b1 = b1f(nsam);
    b2 = b2f(nsam);
    c1 = c1f(a1, b1);
    c2 = c2f(nsam, a1, a2, b2);
    e1 = e1f(a1, c1);
    e2 = e2f(a1, a2, c2);
    return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites -1))) ) ;
    
}

double Correlation(const vector<double>& x, const vector<double>& y)
{
    int n = x.size();
    double ex(0), ey(0), exy(0), ex2(0), ey2(0);
    
    
for (int i = 0; i < n; i++) { //Find the means.
    ex += x[i];
    ey += y[i];
}
ex /= n;
ey /= n;
    double sxx=0; double syy=0; double sxy=0;
for (int i = 0; i < n; i++) { //Compute the correlation coeﬃcient.
    double xt = x[i] - ex;
    double yt = y[i] - ey;
    sxx += xt * xt;
    syy += yt * yt;
    sxy += xt * yt;
}
return sxy/(sqrt(sxx*syy)+1e-20);
}

double littlerSquared(const vector<double>& x, const vector<double>& y)
{
    int n = x.size();
    double ex(0), ey(0), exy(0);
    
    
    for (int i = 0; i < n; i++) { //Find the frequencies.
        ex += x[i];
        ey += y[i];
        exy += (x[i]*y[i]);
    }
    ex /= n;
    ey /= n;
    exy /= n;
 
    double denom = ex * (1-ex)* ey *(1-ey);
    if(denom==0) {return NAN;}
    else{
    double d2 = (exy - ex*ey)*(exy - ex*ey);
    return d2/ denom;
    }
}

double zns(vector< vector< double> > data){
    double sum=0;
    double countCorrs=0;
    for(int i=0; i<data.size();++i){
        for(int j=0;j<data.size();++j){
            if(i<j){
                double r=littlerSquared( data[i], data[j]);
                if(!(r!=r)){
                    sum+=r;
                    countCorrs++;
                }
            }
        }
    }
    
    return sum/countCorrs;
}

double heterozygosity(vector<unsigned> sample){
    
    double pi=0;
    double nd =  sample.size();
    double nnm1 = nd/(nd-1.0) ;
    
    double p1 = (double) accumulate(sample.begin(), sample.end(), 0)/nd;
    pi = 2.0*p1*(1.0 -p1)*nnm1 ;
    
    return(pi);
}

vector<vector<unsigned> > split_by_population(vector<unsigned> sample, vector<unsigned> pop_of_samps){
    vector<vector<unsigned> > splitSNPs;
    
    vector<unsigned> pops = pop_of_samps;
    sort( pops.begin(), pops.end() );
    pops.erase( unique( pops.begin(), pops.end() ), pops.end() );
    unsigned npops = pops.size();
    
    splitSNPs.resize(npops);
    
    for(unsigned i=0; i < sample.size(); ++i){
        unsigned p = pop_of_samps.at(i);
        splitSNPs.at(p).push_back( sample.at(i) );
    }
    return(splitSNPs);
}

vector<double> pi_by_pop(vector<vector<unsigned> > snpvec){
    
    unsigned npops = snpvec.size();
    vector<double> pivec;
    
    for(unsigned p=0; p < npops; ++p){
        double pop_pi = heterozygosity(snpvec.at(p));
        pivec.push_back(pop_pi);
    }
    return(pivec);
}


double calcdxy(vector<vector<unsigned> > samp){
    if(samp.size()!= 2) {return(NAN);}
    
    
    double nx = samp.at(0).size();
    double ny = samp.at(1).size();

    double px = std::accumulate(samp.at(0).begin(), samp.at(0).end(), 0)/nx;
    double py = std::accumulate(samp.at(1).begin(), samp.at(1).end(), 0)/ny;
    
    double correction = (sqrt(nx*ny)-1)/sqrt(nx*ny);
    
    double dxy= correction * (px*(1-py) + py*(1-px));

    return dxy;
}

double fst_nei(double piTot, vector<double> piPop){
    if(piPop.size()!= 2) {return(NAN);}

    if(piTot == 0) return(0);
    double fst = ( piTot - ((piPop.at(0)+piPop.at(1))/2) )/piTot;
    if(fst < 0) return(0);
    return(fst);
}

#endif
