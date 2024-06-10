
//
// Created by Samantha Pope on 3/15/24.
//

#ifndef NEWMSDSCRIPT_VAL_H
#define NEWMSDSCRIPT_VAL_H


#include <cstdio>
#include <string>
#include <sstream>
#include "pointer.h"
#include "Env.h"

using namespace std;
class Expr;
class Env;

CLASS (Val) {
public:
    virtual ~Val() {}
    virtual bool equals (PTR(Val) v)= 0;
    virtual PTR(Val) add_to(PTR(Val) other_val) = 0;
    virtual PTR(Val) mult_with(PTR(Val) other_val) = 0;
    virtual void print(ostream &ostream) = 0;
    virtual bool is_true()=0;
    virtual PTR(Val) call(PTR(Val) actual_arg)=0;

    std::string to_string(){
         std::stringstream st("");
        this->print(st);
        return st.str();
    }
};

class NumVal : public Val{
public:
    int numVal;
    NumVal(int i);
    bool is_true() override;
    bool equals (PTR(Val) v) override;
    PTR(Val) add_to(PTR(Val) other_val) override;
    PTR(Val) mult_with(PTR(Val) other_val) override;
    void print (ostream &ostream) override;
    PTR(Val) call(PTR(Val) actual_arg) override;
};

class BoolVal : public Val{
public:
    bool value;
    bool is_true() override;
    BoolVal(bool passedBool);
    bool equals (PTR(Val) v) override;
    PTR(Val) add_to(PTR(Val) other_val) override;
    PTR(Val) mult_with(PTR(Val) other_val) override;
    void print (ostream &ostream) override;
    PTR(Val) call(PTR(Val) actual_arg) override;

};

class FunVal : public Val {
public:
    string formalArg;
    PTR(Expr) body;
    PTR(Env) env;
    FunVal(std::string formalArgPassed, PTR(Expr) bodyPassed, PTR(Env) env);
    bool equals (PTR(Val) v) override;
    PTR(Val) add_to(PTR(Val) other_val) override;
    PTR(Val) mult_with(PTR(Val) other_val) override;
    void print(std::ostream &ostream) override;
    PTR(Val) call(PTR(Val) actual_arg) override;
    bool is_true() override;
};


#endif //NEWMSDSCRIPT_VAL_H