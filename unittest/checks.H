#ifndef checks_H                                                                
#define checks_H

#include "catch.hpp"
#include "tensor.H"
#include "vector.H"

static Approx approx = Approx::custom().epsilon(0.001);

void check(Foam::vector actual, Foam::vector expected);

void check(Foam::tensor actual, Foam::tensor expected);

#endif
