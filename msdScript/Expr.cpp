#include <sstream>
#include "Expr.h"
#include "val.h"
#include "pointer.h"
#include "Env.h"


/**
 * Constructs a Num object with a given integer value.
 * @param val The integer value to initialize the Num object with.
 */
Num::Num(int val) {
    this->val = val;
}

/**
 * Constructs an Add object with left and right expressions.
 * @param lhs Pointer to the left-hand side expression.
 * @param rhs Pointer to the right-hand side expression.
 */
Add::Add(PTR(Expr) lhs, PTR(Expr) rhs) {
    this->lhs = lhs;
    this->rhs = rhs;
}

/**
 * Constructs a Mult object with left and right expressions.
 * @param lhs Pointer to the left-hand side expression.
 * @param rhs Pointer to the right-hand side expression.
 */
Mult::Mult(PTR(Expr) lhs, PTR(Expr) rhs) {
    this->lhs = lhs;
    this->rhs = rhs;
}

/**
 * Constructs a Var object with a given variable name.
 * @param varPassed The name of the variable.
 */
Var::Var(std::string varPassed) : var(std::move(varPassed)) {
}

/**
 * Constructor for the _Let class
 * @param varName the name of the expression ex: "x"
 * @param expr1 the expression being passed ex: "Num(5)" or "Add(2,3)"
 * @param expr2 the expression being passed ex: "Num(5)" or "Add(2,3)"
 */
_Let::_Let(std::string varName, PTR(Expr) expr1, PTR(Expr) expr2){
    this->varName = varName;
    this->head=expr1;
    this->body=expr2;
}

/**
 * Checks if this Num object is equal to another expression.
 * @param e Pointer to the expression to compare with.
 * @return true if the expressions are equivalent, false otherwise.
 */
bool Num::equals(PTR(Expr) e) {
    PTR(Num) numPtr = CAST(Num)(e);
    return numPtr && this->val == numPtr->val;
}

/**
 * Checks if this Add object is equal to another expression.
 * @param e Pointer to the expression to compare with.
 * @return true if the expressions are equivalent, false otherwise.
 */
bool Add::equals(PTR(Expr) e) {
    PTR(Add) addPtr = CAST(Add)(e);
    return addPtr && ((this->lhs->equals(addPtr->lhs) && this->rhs->equals(addPtr->rhs)) ||
                      (this->lhs->equals(addPtr->rhs) && this->rhs->equals(addPtr->lhs)));
}

/**
 * Checks if this Mult object is equal to another expression.
 * @param e Pointer to the expression to compare with.
 * @return true if the expressions are equivalent, false otherwise.
 */
bool Mult::equals(PTR(Expr) e) {
    PTR(Mult) multPtr = CAST(Mult)(e);
    return multPtr && ((this->lhs->equals(multPtr->lhs) && this->rhs->equals(multPtr->rhs)) ||
                       (this->lhs->equals(multPtr->rhs) && this->rhs->equals(multPtr->lhs)));
}

/**
 * Checks if this Var object is equal to another expression.
 * @param e Pointer to the expression to compare with.
 * @return true if the expressions are equivalent, false otherwise.
 */
bool Var::equals(PTR(Expr) e) {
    PTR(Var) varPtr = CAST(Var)(e);
    return varPtr && this->var == varPtr->var;
}

/**
 * returns if the object passed was the same type of object and if they are the same
 * @param e
 * @return
 */
bool _Let::equals(PTR(Expr) e) {
    auto other = CAST(_Let)(e); //dyamic cast it to that type
    if (other == nullptr) return false; // if type is not the same, return false
    return varName == other->varName && head->equals(other->head) && body->equals(other->body);
}


/**
 * Evaluates the numeric expression.
 * @return The integer value of the Num object.
 */
PTR(Val) Num::interp(PTR(Env) env) {
    return NEW(NumVal)(val);
}

/**
 * Evaluates the addition expression.
 * @return The sum of the left and right expressions.
 */
PTR(Val) Add::interp(PTR(Env) env) {
    return this->lhs->interp(env)->add_to(this->rhs->interp(env));}

/**
 * Evaluates the multiplication expression.
 * @return The product of the left and right expressions.
 */
PTR(Val) Mult::interp(PTR(Env) env) {
    return this->lhs->interp(env)->mult_with(this->rhs->interp(env));
}

/**
 * Evaluates the variable expression. Throws an error because variables cannot be directly evaluated.
 * @throw std::runtime_error when trying to evaluate a variable.
 */
PTR(Val) Var::interp(PTR(Env) env) {
    throw std::runtime_error("no value for variable");
    return NEW(NumVal)(-1);
}

/**
 * Interprets the function by passing in the value the variable is set to and solving
 * @return int that the solution equals
 */
PTR(Val) _Let::interp(PTR(Env) env) {
    PTR(Val) rhsValue = head->interp(env);
    PTR(Env) newEnv = NEW(ExtendedEnv)(varName, rhsValue, env); //TODO is this right? slide 42
    return body->interp(newEnv);
}


/**
 * Prints the Num expression to a given output stream.
 * @param stream The output stream to print to.
 */
void Num::print(std::ostream &stream) {
    stream << std::to_string(val);
}

/**
 * Prints the Add expression to a given output stream.
 * @param stream The output stream to print to.
 */
void Add::print(std::ostream &stream) {
    stream << "(";
    lhs->print(stream);
    stream << "+";
    rhs->print(stream);
    stream << ")";
}

/**
 * Prints the Mult expression to a given output stream.
 * @param stream The output stream to print to.
 */
void Mult::print(std::ostream &stream) {
    stream << "(";
    lhs->print(stream);
    stream << "*";
    rhs->print(stream);
    stream << ")";
}

/**
 * Prints the Var expression to a given output stream.
 * @param stream The output stream to print to.
 */
void Var::print(std::ostream &stream) {
    stream << var;
}

/**
 * Prints out object using a stream
 * @param stream
 */
void _Let::print(std::ostream &stream) {
    stream << "(_let " << varName << "=" << head->to_string() << " _in " << body->to_string() << ")";
}

/**
 * Pretty prints the Num expression with appropriate precedence.
 * @param ot The output stream to print to.
 * @param prec The precedence context in which this expression is being printed.
 */
void Num::pretty_print(std::ostream &ot, precedence_t prec, std::streampos& lastNewLinePos, bool paren) {
    ot << val;
}

/**
 * Pretty prints the Add expression with appropriate precedence.
 * @param ot The output stream to print to.
 * @param prec The precedence context in which this expression is being printed.
 */
void Add::pretty_print(std::ostream &ot, precedence_t prec, std::streampos& lastNewLinePos, bool paren) {
    bool needParens = prec > prec_add;

    //add should always be true??
    if (needParens) ot << "(";
    lhs->pretty_print(ot, static_cast<precedence_t>(prec_add + 1),lastNewLinePos, false);
    ot << " + ";
    rhs->pretty_print(ot, prec_add,lastNewLinePos, needParens);
    if (needParens) ot << ")";
}

/**
 * Pretty prints the Mult expression with appropriate precedence.
 * @param ot The output stream to print to.
 * @param prec The precedence context in which this expression is being printed.
 */
void Mult::pretty_print(std::ostream &ot, precedence_t prec, std::streampos& lastNewLinePos, bool paren) {
    bool needParens = prec > prec_mult;
    if (needParens) ot << "(";
    this->lhs->pretty_print(ot, static_cast<precedence_t>(prec_mult + 1),lastNewLinePos, needParens);
    ot << " * ";
    this->rhs->pretty_print(ot, prec_mult,lastNewLinePos, needParens);
    if (needParens) ot << ")";
}

/**
 * Pretty prints the Var expression with appropriate precedence.
 * @param ot The output stream to print to.
 * @param prec The precedence context in which this expression is being printed.
 */
void Var::pretty_print(std::ostream &ot, precedence_t prec, std::streampos& lastNewLinePos, bool paren) {
    ot << var;
}


/**
 * Prints out object in a more readable way using precedence, stream and streampos
 * @param ot
 * @param prec
 * @param lastNewLinePos
 */
void _Let::pretty_print(std::ostream &ot, precedence_t prec, std::streampos& lastNewLinePos, bool paren) {
    if (!paren && prec != prec_none) {
        ot << "(";
    }
    std::streampos letPosition = ot.tellp();
    std::streampos depth = letPosition - lastNewLinePos;
    ot << "_let " << this->varName<<" = ";
    //print bound expression with passing difference
    this->head->pretty_print(ot,prec_none, depth,paren);
    ot << "\n ";
    std::streampos nextPos = ot.tellp();
    //start print with indentation
    for ( int i = 0; i < letPosition - lastNewLinePos; i++ ) {
        ot << " ";
    }
    ot<< "_in  ";
    this->body->pretty_print(ot, prec_none, nextPos,paren);
    if (!paren && prec != prec_none) {
        ot << ")";
    }

}

void Num::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp(); //initiate to 0
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}

void Var::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp(); //initiate to 0
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);

}

void Mult::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos = ot.tellp(); //initiate to 0
    this-> pretty_print(ot,prec_mult,lastNewLinePos, false);
}

void Add::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_add,lastNewLinePos, false);
}

void _Let::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}


 //-----------------------EqExpr--------------------------
/**
 * Constructor for EqExpr
 * @param lhs
 * @param rhs
 */
EqExpr::EqExpr(PTR(Expr) lhs, PTR(Expr) rhs) {
    this->rhs = rhs;
    this->lhs =lhs;
}

bool EqExpr::equals(PTR(Expr) e) {
    PTR(EqExpr) eq = CAST(EqExpr)(e);
    if (eq == nullptr) return false;
    return lhs->equals(eq->lhs) && rhs->equals(eq->rhs);
}

PTR(Val) EqExpr::interp(PTR(Env) env) {
    return NEW(BoolVal)(rhs->interp(env)->equals(lhs->interp(env)));
}



void EqExpr::print(std::ostream &stream) {
    lhs->print(stream);
    stream << " == ";
    rhs->print(stream);
}

void EqExpr::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_add,lastNewLinePos, false);
}

void EqExpr::pretty_print(std::ostream &ot, precedence_t prec, std::streampos &lastNewLinePos, bool paren) {
    if (paren) ot << "(";
    lhs->pretty_print(ot, prec_add, lastNewLinePos, false);
    ot << " == ";
    rhs->pretty_print(ot, prec_add, lastNewLinePos, false);
    if (paren) ot << ")";
}


//----------------------IfExpr-----------------------

IfExpr::IfExpr(PTR(Expr) condition, PTR(Expr) thenExpr, PTR(Expr) elseExpr) {
    this->condition=condition;
    this->thenExpr=thenExpr;
    this->elseExpr=elseExpr;
}

bool IfExpr::equals(PTR(Expr) e) {
    PTR(IfExpr) other = CAST(IfExpr)(e);
    return other != nullptr
           && condition->equals(other->condition)
           && thenExpr->equals(other->thenExpr)
           && elseExpr->equals(other->elseExpr);
}

PTR(Val) IfExpr::interp(PTR(Env) env) {
    PTR(Val) condVal = condition->interp(env);
    PTR(BoolVal) boolVal = CAST(BoolVal)(condVal);
    if (boolVal == nullptr) {
        throw std::runtime_error("Condition expression did not evaluate to a boolean");
    }
    PTR(Val) result = boolVal->is_true() ? thenExpr->interp(env) : elseExpr->interp(env);
    return result;
}




void IfExpr::print(std::ostream &stream) {
    stream << "_if ";
    condition->print(stream);
    stream << " _then ";
    thenExpr->print(stream);
    stream << " _else ";
    elseExpr->print(stream);
}

void IfExpr::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}

void IfExpr::pretty_print(std::ostream &ot, precedence_t prec, std::streampos &lastNewLinePos, bool paren) {
    ot << "_if ";
    condition->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
    lastNewLinePos = ot.tellp();
    ot << "_then ";
    thenExpr->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
    ot << "_else ";
    lastNewLinePos = ot.tellp();
    elseExpr->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
}

//------------BoolExpr---------------------

BoolExpr::BoolExpr(bool value) {
    this->value =value;
}

bool BoolExpr::equals(PTR(Expr) e) {
    PTR(BoolExpr) be = CAST(BoolExpr)(e);
    return be != nullptr && value == be->value;
}

PTR(Val) BoolExpr::interp(PTR(Env) env) {
    return NEW(BoolVal)(value);
}



void BoolExpr::print(std::ostream &stream) {
    stream << (value ? "_true" : "_false");
}

void BoolExpr::pretty_print(std::ostream &ot, precedence_t prec, std::streampos &lastNewLinePos, bool paren) {
    print(ot);
}

void BoolExpr::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}

FunExpr::FunExpr(std::string passedArg, PTR(Expr) passedBody) {
    this->formalArg = passedArg;
    this->body = passedBody;
}


bool FunExpr::equals(PTR(Expr) e) {
    PTR(FunExpr) funPtr = CAST(FunExpr)(e);
    if (funPtr == nullptr){
        return false;
    }
    return this->formalArg == funPtr->formalArg && this->body->equals(funPtr->body);
}

PTR(Val) FunExpr::interp(PTR(Env) env) {
    if (env == nullptr){
        env = Env::empty;
    }
    return NEW(FunVal)(formalArg, body, env);
}


void FunExpr::print(ostream& o){
    o << "(_fun (" << this->formalArg << ") " << this->body->to_string() << ")";
}

void FunExpr::pretty_print(std::ostream &ot, precedence_t prec, std::streampos &lastNewLinePos, bool paren) { //OPTIONAL so i did it how i want
    ot << "_fun ";
    ot << formalArg;
    lastNewLinePos = ot.tellp();
    ot << "==";
    body->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
}

void FunExpr::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}

CallExpr::CallExpr(PTR(Expr) toBeCalled, PTR(Expr) actualArg){
    this->toBeCalled = toBeCalled;
    this->actualArg = actualArg;
};

bool CallExpr::equals(PTR(Expr) e){
    PTR(CallExpr) callPtr = CAST(CallExpr)(e);
    if (callPtr == nullptr){
        return false;
    }
    return this->toBeCalled->equals(callPtr->toBeCalled) && this->actualArg->equals(callPtr->actualArg);
}
PTR(Val) CallExpr::interp(PTR(Env) env){
    return this->toBeCalled->interp(env)->call(actualArg->interp(env));
}

void CallExpr::print(ostream &ostream){
    ostream << "(" << this->toBeCalled->to_string() << ") (" << this->actualArg->to_string() << ")";
}

void CallExpr::pretty_print(std::ostream &ot, precedence_t prec, std::streampos &lastNewLinePos, bool paren) { //optional
    ot << "to be call: ";
    toBeCalled->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
    lastNewLinePos = ot.tellp();
    ot << "actual arg: ";
    actualArg->pretty_print(ot, prec_none, lastNewLinePos, false);
    ot << "\n";
}

void CallExpr::pretty_print_at(std::ostream &ot) {
    std::streampos lastNewLinePos =ot.tellp();
    this-> pretty_print(ot,prec_none,lastNewLinePos, false);
}