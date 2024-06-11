//
// Created by Samantha Pope on 4/3/24.
//

#include "Env.h"

PTR(Env) Env::empty = NEW(EmptyEnv)();

PTR(Val) EmptyEnv::lookup(std::string findName) {
    throw std::runtime_error("free variable: " + findName);
};

ExtendedEnv::ExtendedEnv(std::string name_, PTR(Val) val_, PTR(Env) rest_) {
    name = name_;
    val = val_;
    rest = rest_;
}

PTR(Val) ExtendedEnv::lookup(std::string findName) {
    if(findName == name){
        return val;
    } else {
        return rest->lookup(findName);
    }
};