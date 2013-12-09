#ifndef _MPS_OCTAVE_SUPPORT_H__
#define _MPS_OCTAVE_SUPPORT_H__

#include <set>
#include <string>

namespace mpsolve {
  namespace octave {
    
    class Algorithms {
      
    public:
      
      /**
       * @brief Return true if the algorithm specified by the given string is supported. 
       */ 
      static bool isSupported (std::string alg);
      
    private:
      
      static std::set<std::string> m_supportedAlgs; 
      static bool m_initPerformed; 
      
      static void initSupportedAlgs (); 
      
    }; 
  
  } 
}

#endif
