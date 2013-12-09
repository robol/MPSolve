#include "octave_support.h"

using namespace mpsolve::octave; 

std::set<std::string> Algorithms::m_supportedAlgs;
bool Algorithms::m_initPerformed = false; 

void
Algorithms::initSupportedAlgs () {
  if (! m_initPerformed) {
    m_initPerformed = true; 
    m_supportedAlgs.insert (std::string("s")); 
    m_supportedAlgs.insert (std::string("h")); 
  }
}

bool
Algorithms::isSupported (std::string alg) {
  initSupportedAlgs (); 
  return m_supportedAlgs.find(alg) != m_supportedAlgs.end(); 
}
