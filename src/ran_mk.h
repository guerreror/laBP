/*
 * ran_mk
 *  InvCoal
 *  Version:  labp_v17
 *
 * Routines for generating various kinds of random numbers.
 *
 * Started by Shane Pope & Mark Kirkpatrick, August 2008
 *
 * Modified by Rafael Guerrero 2009-2016
 * Included exponential distributed numbers, switched references to Boost -> std library
 *
 */

#ifndef RANDOM_
#define RANDOM_


#include <random>

#include <boost/random.hpp>
using boost::variate_generator;

extern std::mt19937_64 gen;	// Declared externally in main.cpp so that it can be initialized there and called in random.cpp

typedef variate_generator<std::mt19937&, std::uniform_int_distribution<> > rng_uniform_int_t;
typedef variate_generator<std::mt19937&, std::uniform_real_distribution<> > rng_uniform_real_t;
typedef variate_generator<std::mt19937&, std::exponential_distribution<> > rng_exponential_t;
typedef variate_generator<std::mt19937&, std::binomial_distribution<> > rng_binomial_t;



template<class PRNG, class Dist>
inline variate_generator<PRNG&, Dist> make_gen(PRNG & rng, Dist d)
{
  return variate_generator<PRNG&, Dist>(rng, d);
}



inline rng_uniform_int_t::result_type randint(int lower_bound, int upper_bound)
{
//
// randint(a, b) returns a uniformly distributed random integer on [a, b]
//
    return make_gen(gen, std::uniform_int_distribution<>(lower_bound, upper_bound))();
}



inline rng_uniform_real_t::result_type randreal(double lower_bound = 0, double upper_bound = 1)
{
//
// Returns a uniformly distributed random deviate.  
//	randreal() is uniform on [0, 1]
//	randreal(a, b) is uniform on [a, b].

    return make_gen(gen, std::uniform_real_distribution<>(lower_bound, upper_bound))();
}



inline rng_exponential_t::result_type randexp(double lambda)
{
//
// Returns a random deviate from an exponential distribution with parameter lambda
//
    return make_gen(gen, std::exponential_distribution<>(lambda))();
}

inline int randpois( double n ) {
    using dist_t = std::poisson_distribution<int>;
    static dist_t d{n};
    return d(gen);
}

inline rng_binomial_t::result_type randbinom(unsigned int trials, double p)
{
    //
    // Returns a random deviate from an poisson distribution with parameter lambda
    //
    return make_gen(gen, std::binomial_distribution<>(trials, p))();
}

#endif //RANDOM_
