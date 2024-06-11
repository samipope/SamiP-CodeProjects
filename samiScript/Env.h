//
// Created by Samantha Pope on 4/3/24.
//

#ifndef NEWMSDSCRIPT_ENV_H
#define NEWMSDSCRIPT_ENV_H
#include "pointer.h"
#include "val.h"
#include "Expr.h"
#include "Env.h"

class Val;
class Expr;

CLASS(Env) {
public:
    static PTR(Env) empty;
    virtual PTR(Val) lookup(std::string findName) =0;
    virtual ~Env() {} //virtual destructor
};

class EmptyEnv : public Env {
public:
    PTR(Val) lookup(std::string findName) override;
};

class ExtendedEnv : public Env {
    std::string name;
    PTR(Val) val;
    PTR(Env) rest;

public:
    ExtendedEnv(std::string name_, PTR(Val) val_, PTR(Env) rest_);
    PTR(Val) lookup(std::string findName) override;
};


#endif //NEWMSDSCRIPT_ENV_H
