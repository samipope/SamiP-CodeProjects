/// Title: Parse.cpp
/// Author: Samantha Pope
/// Date: 02.14.2024
/// Scope: parses user input for our script to understand.
/// Handles all inputs/object types that our user could pass in

#include "Expr.h"
#include "val.h"
#include <iostream>
#include <istream>
#include "parse.h"
#include "pointer.h"
using namespace std;

/**
 * Consumes a specific word from the input stream.
 * If the word does not match exactly, throws a runtime error.
 *
 * @param in The input stream to consume from.
 * @param str The string that needs to be consumed from the input stream.
 * @throws runtime_error If the actual input does not match the expected string.
 */
void consume_word(istream &in, string str){
    for(char c : str){
        if (in.get()!=c){
            throw runtime_error("consume mismatch");
        }
    }
}

/**
 * Parses an if-expression from the input stream.
 * Constructs an IfExpr object based on the conditional, then, and else expressions read from the stream.
 *
 * @param stream The input stream to parse from.
 * @return A pointer to the constructed IfExpr object.
 */
PTR(Expr) parse_if( std::istream &stream ) {
    skip_whitespace(stream);
    PTR(Expr) ifStatement = parse_expr(stream);
    skip_whitespace(stream);
    consume_word(stream, "_then");
    skip_whitespace(stream);
    PTR(Expr) thenStatement = parse_expr(stream);
    skip_whitespace(stream);
    consume_word(stream, "_else");
    skip_whitespace(stream);
    PTR(Expr) elseStatement = parse_expr(stream);
    return NEW(IfExpr)(ifStatement, thenStatement, elseStatement);
}

/**
 * Parses an expression from the input stream.
 * Handles equality and passes control to parse_comparg for arithmetic operations.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed Expr object, representing the parsed expression.
 */
PTR(Expr) parse_expr(std::istream &in) {
    PTR(Expr) e = parse_comparg(in);
    skip_whitespace(in);
    if (in.peek() == '='){
        consume(in, '=');
        if (in.peek() != '='){
            throw runtime_error("need '=='!");
        }
        consume(in, '=');
        PTR(Expr) rhs = parse_expr(in);
        return NEW (EqExpr)(e, rhs);
    }
    return e;
}

/**
 * Parses an addition or subtraction expression from the input stream.
 * Constructs an Add object if an addition operation is detected.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed Expr object, either the parsed expression or an addition expression.
 */
PTR(Expr) parse_comparg(istream &in){
    PTR(Expr) e = parse_addend(in);
    skip_whitespace(in);
    if (in.peek() == '+'){
        consume(in, '+');
        PTR(Expr) rhs = parse_comparg(in);
        return NEW(Add)(e, rhs);
    }
    return e;
}

/**
 * Parses a multiplication expression from the input stream.
 * Constructs a Mult object if a multiplication operation is detected.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed Expr object, either the parsed expression or a multiplication expression.
 */
PTR(Expr) parse_addend(std::istream &in) {
    PTR(Expr) e;
    e = parse_multicand(in);
    skip_whitespace(in);

    int c = in.peek();
    if (c == '*') {
        consume(in, '*');
        skip_whitespace(in);
        PTR(Expr) rhs = parse_addend(in);
        return NEW(Mult)(e, rhs);
    } else {
        return e;
    }
}

/**
 * Parses a term from the input stream. A term is a sequence of letters.
 *
 * @param in The input stream to parse from.
 * @return The parsed term as a string.
 */
string parse_term(istream &in){
    string term;
    while (true) {
        int letter = in.peek();
        if (isalpha(letter)) {
            consume(in, letter);
            char c = letter;
            term += c;
        }
        else
            break;
    }
    return term;
}

/**
 * Parses an expression that can be a function call from the input stream.
 * Handles function call syntax with parentheses.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed Expr object, possibly representing a function call.
 */
PTR(Expr) parse_multicand(istream &in) {
    PTR(Expr) e = parse_inner(in);
    while (in.peek() == '(') {
        consume(in, '(');
        PTR(Expr) actual_arg = parse_expr(in);
        consume(in, ')');
        e = NEW(CallExpr)(e, actual_arg);
    }
    return e;
}

/**
 * Parses the most basic expressions from the input stream, including numbers, variables, and expressions in parentheses.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed Expr object, representing the parsed inner expression.
 */
PTR(Expr) parse_inner(std::istream &in) {
    skip_whitespace(in);
    int c = in.peek();

    if ((c == '-') || isdigit(c)){
        return parse_num(in);
    }

    else if (c == '(') {
        consume(in, '(');
        PTR(Expr) e = parse_comparg(in);
        skip_whitespace(in);
        c = in.get();
        if (c != ')'){
            throw runtime_error("missing closing parentheses");
        }
        return e;
    }

    else if (isalpha(c)) {
        return parse_var(in);
    }

    else if (c=='_'){
        consume(in, '_');

        string term = parse_term(in);

        if(term == "let"){
            return parse_let(in);
        }
        else if(term == "if"){
            return parse_if(in);
        }
        else if(term == "true"){
            return NEW(BoolExpr)(true);
        }
        else if(term == "false"){
            return NEW(BoolExpr)(false);
        }
        else if(term == "fun"){
            return parse_fun(in);
        }
        else{
            throw runtime_error("invalid input");
        }
    }
    else {
        consume(in, c);
        throw runtime_error("invalid input");
    }
}

/**
 * Parses a number from the input stream. Handles negative numbers.
 *
 * @param in The input stream to parse from.
 * @return A pointer to a Num object, representing the parsed number.
 */
PTR(Expr) parse_num(std::istream &in) {

    int n = 0;
    bool negative = false;
    bool digitSeen = false;

    if (in.peek() == '-') {
        negative = true;
        consume(in, '-');
    }

    while (1) {
        int c = in.peek();
        if (isdigit(c)) {
            consume(in, c);
            n = n * 10 + (c - '0');
            digitSeen = true;
        } else
            break;
    }
    if (negative && !digitSeen){
        throw std::runtime_error("invalid input");
    }
    if (negative) {
        n = -n;
    }
    return NEW(Num)(n);
}
/**
 * Consumes a specific character from the input stream. If the character does not match, throws a runtime error.
 *
 * @param in The input stream to consume from.
 * @param expect The character expected to be consumed.
 * @throws runtime_error If the actual character does not match the expected character.
 */
void consume(std::istream &in, int expect) {
    int c = in.get();
    if (c != expect) {
        throw std::runtime_error("consume mismatch");
    }
}

/**
 * Consumes a specific string from the input stream. If the string does not match exactly, throws a runtime error.
 *
 * @param stream The input stream to consume from.
 * @param str The string that needs to be consumed from the input stream.
 * @throws runtime_error If the actual input does not match the expected string.
 */
void consume( std::istream & stream, const std::string & str)
{
    for ( char expect : str )
    {
        const int c = stream.get();
        if ( c != expect )
            throw std::runtime_error( "consume(): mismatch" );
    }
}

/**
 * Skips whitespace characters in the input stream.
 *
 * @param in The input stream to skip whitespace in.
 */
void skip_whitespace(std::istream &in) {
    while (1) {
        int c = in.peek();
        if (!isspace(c))
            break;
        consume(in, c);
    }
}

/**
 * Parses an entire expression from the input stream.
 * This function ensures that the entire input is consumed by the end of the parsing process.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the parsed Expr object.
 * @throws runtime_error If there is any remaining input after parsing the expression.
 */
PTR(Expr) parse(std::istream &in) {
    PTR(Expr) e;
    e = parse_expr(in);
    skip_whitespace(in);
    if (!in.eof()) {
        throw std::runtime_error("invalid input");
    }
    return e;
}

/**
 * Parses an input expression from standard input and constructs its corresponding expression tree.
 * @return A pointer to the constructed Expr object, representing the parsed expression.
 */
PTR(Expr) parseInput() {
    std::string input;
    getline(std::cin, input);
    std::cout << "input : " << input << std::endl;
    std::stringstream ss(input);
    return parse_comparg(ss);
}

/**
 * Parses a let-expression from the input stream.
 * Constructs a _Let object based on the parsed variable, value, and body expressions.
 *
 * @param in The input stream to parse from.
 * @return A pointer to the constructed _Let object.
 */
PTR(Expr) parse_let(std::istream &in){
    skip_whitespace(in);
    PTR(Expr) e = parse_var(in);
    string lhs = e->to_string();
    skip_whitespace(in);
    consume(in, '=');
    skip_whitespace(in);
    PTR(Expr) rhs = parse_comparg(in);
    skip_whitespace(in);
    consume_word(in, "_in");
    skip_whitespace(in);
    PTR(Expr) body = parse_comparg(in);
    return NEW(_Let)(lhs, rhs, body);
}


/**
 * Parses a variable name from the input stream and creates a variable expression.
 * The function reads characters from the input stream until it encounters a character that is not a letter.
 * It then constructs a Var object representing a variable expression with the parsed name.
 *
 * @param in The input stream to parse from.
 * @return A pointer to a Var object, representing the parsed variable expression.
 */
PTR(Expr)parse_var(std::istream &in){
    std::string var;
    while(true){
        int c = in.peek();
        if (isalpha(c)){
            consume(in, c);
            var += static_cast<char>(c);
        } else {
            break;
        }
    }
    return NEW(Var)(var);
}

/**
 * Parses an expression from a string.
 * This function wraps the given string into an input stream and delegates the parsing to the `parse` function,
 * which constructs an expression tree based on the input.
 *
 * @param s The string containing the expression to parse.
 * @return A pointer to an Expr object, representing the root of the parsed expression tree.
 */
PTR(Expr) parse_str(const string& s){
    istringstream in(s);
    return parse (in);
}

/**
 * Parses a function expression from the input stream.
 * This function expects the input stream to contain a function declaration, starting with an opening parenthesis,
 * followed by a variable name (the function parameter), a closing parenthesis, and the function body as an expression.
 * It constructs a FunExpr object representing the parsed function expression.
 *
 * @param in The input stream to parse from.
 * @return A pointer to a FunExpr object, representing the parsed function expression.
 */
PTR(Expr) parse_fun(istream &in){
    skip_whitespace(in);
    consume(in, '(');
    PTR(Expr) e = parse_var(in);
    string var = e->to_string();
    consume(in, ')');
    skip_whitespace(in);
    e = parse_expr(in);
    return NEW(FunExpr)(var, e);
}
