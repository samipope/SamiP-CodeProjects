//
// Created by Samantha Pope on 3/15/24.
//

#include "val.h"

#include <utility>
#include "Expr.h"



NumVal::NumVal(int i) {
    numVal = i;
}


bool NumVal::equals(PTR(Val) v) {
    //Insert implementation
    PTR(NumVal) numPointer = CAST(NumVal)(v);
    if (numPointer == nullptr){
        return false;
    }
    return this->numVal == numPointer ->numVal;
}

PTR(Val) NumVal::add_to(PTR(Val) other_val) {
    PTR(NumVal) other_num = CAST(NumVal)(other_val);
    if (other_num == nullptr) throw runtime_error("You can't add a non-number!");
    return NEW(NumVal)( (unsigned) other_num->numVal + (unsigned) this->numVal); //assigned to unsigned help prevent undefined behavior
}

PTR(Val) NumVal::mult_with(PTR(Val) other_val) {
    PTR(NumVal) other_num = CAST(NumVal)(other_val);
    if(other_num == nullptr) throw runtime_error("You can't mult a non-number!");
    return NEW(NumVal)((unsigned)this->numVal * (unsigned)other_num->numVal); //assigned to unsigned help prevent undefined behavior
}

void NumVal::print(std::ostream &ostream) {
    ostream << numVal;
}

bool NumVal::is_true() {
    throw std::runtime_error("cannot use is_true on NumVal");
}

PTR(Val) NumVal::call(PTR(Val) actualArg){
    throw runtime_error("Cannot call NumVal!");
}


BoolVal::BoolVal(bool passedBool){
    value = passedBool;
}

bool BoolVal::equals(PTR(Val) v) {
    PTR(BoolVal) bv = CAST(BoolVal)(v);
    return bv != nullptr && value == bv->value;
}



PTR(Val) BoolVal::add_to(PTR(Val) other_val) {
    throw std::runtime_error("Cannot add boolean values");
}

PTR(Val) BoolVal::mult_with(PTR(Val) other_val) {
    throw std::runtime_error("Cannot multiply boolean values");
}

void BoolVal::print(std::ostream &ostream) {
    ostream <<:: to_string(value);
}

bool BoolVal::is_true() {
    return value;
}

PTR(Val) BoolVal::call(PTR(Val) actualArg){
    throw runtime_error("Cannot call BoolVal");
}

FunVal::FunVal(std::string formalArgPassed, PTR(Expr) bodyPassed, PTR(Env) env) {
    if(env == nullptr){
        env = Env::empty;
    }
    this->formalArg = std::move(formalArgPassed);
    this->body = std::move(bodyPassed);
}


bool FunVal::equals(PTR(Val) v) {
    PTR(FunVal) funPtr = CAST(FunVal)(v);
    if (funPtr == nullptr){
        return false;
    }
    return this->formalArg == funPtr->formalArg && this->body->equals(funPtr->body);
}

PTR(Val) FunVal::add_to(PTR(Val) other_val) {
    throw runtime_error("can't add");
}

PTR(Val) FunVal::mult_with(PTR(Val) other_val) {
    throw runtime_error("can't multiply");
}

void FunVal::print(std::ostream &ostream) {
}

bool FunVal::is_true() {
    return false;
}

PTR(Val) FunVal::call(PTR(Val) actual_arg) {
    return body->interp(NEW(ExtendedEnv)(formalArg,actual_arg,env));
}