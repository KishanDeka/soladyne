#include "solar_dynamo/model.hpp"
#include <cassert>
#include <cmath>
int main(){solar_dynamo::Parameters p;p.relaxed_initial_state=false;solar_dynamo::Model m(p);m.initialize("ignored.dat");auto s=m.snapshot();assert(s.omega.size()==257u*257u);m.step(2);assert(m.step_number()==2);assert(std::abs(m.time()-2*p.dt)<1e-15);}
