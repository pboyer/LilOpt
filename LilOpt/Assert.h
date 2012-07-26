//
//  Error.h
//  LilOpt
//
//  Created by Peter Boyer on 7/25/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilLM_Assert_h
#define LilLM_Assert_h

#include <iostream>

namespace LilOpt {
    
    void AssertOrExit(bool predicate, std::string message) {
    
        if (!predicate) {
        
            std::cerr << "LILOPT RUNTIME ERROR" << std::endl;
            std::cerr << "====================" << std::endl;
            std::cerr << message << std::endl;
            std::cerr << "EXITING..." << std::endl;
            exit(1);
        
        }
    
    }
}

#endif
