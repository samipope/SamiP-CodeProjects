/// Title: ExprTest.cpp
/// Author: Samantha Pope
/// Date: 1.22.2024
/// Scope: This document includes all of my tests for my script using the Catch2 testing format.
///


#include "catch.h"
#include "Expr.h"
#include "parse.h"
#include "cmdline.h"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include "val.h"
#include "pointer.h"

TEST_CASE("Num equals tests", "[Num]") {
    auto num1 = NEW(Num)(5);
    auto num2 = NEW(Num)(5);
    auto num3 = NEW(Num)(10);
    SECTION("Equal numbers") {
        REQUIRE(num1->equals(num2));
    }
    SECTION("Not equal numbers") {
        REQUIRE_FALSE(num1->equals(num3));
    }
    SECTION("Comparison with null") {
        REQUIRE_FALSE(num1->equals(nullptr));
    }
}


TEST_CASE("Add equals tests", "Add") {
    auto n1 = NEW(Num)(2);
    auto n2 = NEW(Num)(5);
    auto n3 = NEW(Num)(2);
    auto n4 = NEW(Num)(5);
    auto mult1 = NEW(Mult)(n1, n2);
    auto mult2 = NEW(Mult)(n3, n4);
    auto mult3 = NEW(Mult)(n2, n1);  // Different order
    auto differentNum = NEW(Num)(9);
    auto differentMult = NEW(Mult)(n1, differentNum);
    SECTION("Equal Mults") {
        REQUIRE(mult1->equals(mult2));
    }
    SECTION("Equal Mults with different order") {
        REQUIRE(mult1->equals(mult3));
    }
    SECTION("Not equal Mults") {
        REQUIRE_FALSE(mult1->equals(differentMult));
    }
}

TEST_CASE("Var equals tests", "[Var]") {
    auto varExpr1 = NEW(Var)("x");
    auto varExpr2 = NEW(Var)("x");
    auto varExpr3 = NEW(Var)("y");
    auto numExpr = NEW(Num)(5);
    SECTION("Equal VarExprs with same name") {
        REQUIRE(varExpr1->equals(varExpr2));
    }
    SECTION("Not Equal VarExprs with different names") {
        REQUIRE_FALSE(varExpr1->equals(varExpr3));
    }
    SECTION("Not Equal when compared with different type (Num)") {
        REQUIRE_FALSE(varExpr1->equals(numExpr));
    }
}

TEST_CASE("interp tests", "All Expressions") {
    SECTION("Num interp") {
        CHECK(NEW(Num)(3)->interp(Env::empty)->equals(NEW(NumVal)(3)));
        CHECK(NEW(Num)(5)->interp(Env::empty)->equals(NEW(NumVal)(5)));
        CHECK(NEW(Num)(-18)->interp(Env::empty)->equals(NEW(NumVal)(-18)));
        CHECK(NEW(Num)(-3)->interp(Env::empty)->equals(NEW(NumVal)(-3)));
        CHECK(NEW(Num)(0)->interp(Env::empty)->equals(NEW(NumVal)(0)));
    }

    SECTION("Add interp") {
        CHECK(NEW(Add)(NEW(Num)(3), NEW(Num)(2))->interp(Env::empty)->equals(NEW(NumVal)(5)));
        CHECK(NEW(Add)(NEW(Num)(5), NEW(Num)(-4))->interp(Env::empty)->equals(NEW(NumVal)(1)));
        CHECK(NEW(Add)(NEW(Num)(-3), NEW(Num)(3))->interp(Env::empty)->equals(NEW(NumVal)(0)));
        CHECK(NEW(Add)(NEW(Num)(-3), NEW(Num)(-3))->interp(Env::empty)->equals(NEW(NumVal)(-6)));
        CHECK(NEW(Add)(NEW(Num)(0), NEW(Num)(10))->interp(Env::empty)->equals(NEW(NumVal)(10)));
    }

    SECTION("Mult interp") {
        CHECK(NEW(Mult)(NEW(Num)(3), NEW(Num)(2))->interp(Env::empty)->equals(NEW(NumVal)(6)));
        CHECK(NEW(Mult)(NEW(Num)(5), NEW(Num)(4))->interp(Env::empty)->equals(NEW(NumVal)(20)));
        CHECK(NEW(Mult)(NEW(Num)(-3), NEW(Num)(6))->interp(Env::empty)->equals(NEW(NumVal)(-18)));
        CHECK(NEW(Mult)(NEW(Num)(-3), NEW(Num)(-3))->interp(Env::empty)->equals(NEW(NumVal)(9)));
        CHECK(NEW(Mult)(NEW(Num)(0), NEW(Num)(10))->interp(Env::empty)->equals(NEW(NumVal)(0)));
    }


}
TEST_CASE("subst tests", "All Expressions") {
auto num5 = NEW(Num)(5);
auto num10 = NEW(Num)(10);
auto varX = NEW(Var)("x");
auto varY = NEW(Var)("y");

}

TEST_CASE("to_string tests", "all expressions") {
    CHECK(NEW(Var)("x")->to_string() == "x");
    CHECK(NEW(Add)(NEW(Num)(1), NEW(Num)(2))->to_string() == "(1+2)");
    CHECK(NEW(Mult)(NEW(Num)(3), NEW(Num)(4))->to_string() == "(3*4)");
    CHECK(NEW(Num)(-5)->to_string() == "-5");
    CHECK(NEW(Num)(0)->to_string() == "0");
    CHECK(NEW(Add)(NEW(Num)(1), NEW(Num)(-2))->to_string() == "(1+-2)");
    CHECK(NEW(Mult)(NEW(Num)(0), NEW(Num)(4))->to_string() == "(0*4)");
    CHECK(NEW(Mult)(NEW(Num)(3), NEW(Num)(-4))->to_string() == "(3*-4)");
    CHECK(NEW(Add)(NEW(Num)(-1), NEW(Mult)(NEW(Num)(2), NEW(Num)(-3)))->to_string() == "(-1+(2*-3))");
    CHECK(NEW(Add)(NEW(Var)("x"), NEW(Mult)(NEW(Add)(NEW(Num)(0), NEW(Num)(-5)), NEW(Var)("y")))->to_string() == "(x+((0+-5)*y))");
}



TEST_CASE("pretty_print_at Tests", "All expression classes") {
    std::stringstream ss3;
    NEW(Add)(NEW(Num)(4), NEW(Num)(5))->pretty_print_at(ss3);
    CHECK(ss3.str() == "4 + 5");

    std::stringstream ss4;
    NEW(Add)(NEW(Num)(4), NEW(Mult)(NEW(Num)(5), NEW(Num)(6)))->pretty_print_at(ss4);
    CHECK(ss4.str() == "4 + 5 * 6");

    std::stringstream ss5;
    NEW(Add)(NEW(Mult)(NEW(Num)(1), NEW(Num)(2)), NEW(Add)(NEW(Num)(3), NEW(Mult)(NEW(Num)(4), NEW(Num)(5))))->pretty_print_at(ss5);
    CHECK(ss5.str() == "1 * 2 + 3 + 4 * 5");
}

TEST_CASE("_Let Tests"){

    SECTION("equals"){
        auto let1 = NEW(_Let)("x", NEW(Num)(5), NEW(Var)("x"));
        auto let2 = NEW(_Let)("x", NEW(Num)(5), NEW(Var)("x"));
        auto letDiffVar = NEW(_Let)("y", NEW(Num)(5), NEW(Var)("x"));
        auto letDiffExpr = NEW(_Let)("x", NEW(Num)(5), NEW(Var)("y"));
        auto letCommutative = NEW(_Let)("x", NEW(Mult)(NEW(Num)(2), NEW(Num)(3)), NEW(Var)("x"));
        auto letCommutativeDiffOrder = NEW(_Let)("x", NEW(Mult)(NEW(Num)(3), NEW(Num)(2)), NEW(Var)("x"));

        CHECK(let1->equals(let2));
        CHECK_FALSE(let1->equals(letDiffVar));
        CHECK(let1->equals(letDiffExpr) == false);
        CHECK(letCommutative->equals(letCommutativeDiffOrder));
    }

}
TEST_CASE("parse") {
    SECTION("parsing single numbers") {
        CHECK(parse_str("1")->equals(NEW(Num)(1)));
        CHECK(parse_str("(((1)))")->equals(NEW(Num)(1)));
        CHECK(parse_str("  \n 5  ")->equals(NEW(Num)(5)));
    }

    SECTION("Handling invalid input") {
        CHECK_THROWS_WITH(parse_str("(-5"), "missing closing parentheses");
        CHECK_THROWS_WITH(parse_str("-"), "invalid input");
        CHECK_THROWS_WITH(parse_str(" -   "), "invalid input");
        CHECK_THROWS_WITH(parse_str("x_z"), "invalid input");
    }

    SECTION("Parsing variables") {
        CHECK(parse_str("xyz")->equals(NEW(Var)("xyz")));
        CHECK(parse_str("xYz")->equals(NEW(Var)("xYz")));
    }

    SECTION("Parsing addition expressions") {
        CHECK(parse_str("x + y")->equals(NEW(Add)(NEW(Var)("x"), NEW(Var)("y"))));
    }

    SECTION("Parsing multiplication expressions") {
        CHECK(parse_str("x * y")->equals(NEW(Mult)(NEW(Var)("x"), NEW(Var)("y"))));
    }

    SECTION("Parsing negative numbers") {
        CHECK(parse_str("-5")->equals(NEW(Num)(-5)));
        CHECK_THROWS_WITH(parse_str(" -   5"), "invalid input");
    }
}

TEST_CASE("NumVal") {
    SECTION("equals method") {
        auto a = NEW(NumVal)(10);
        auto b = NEW(NumVal)(10);
        auto c = NEW(NumVal)(5);
        CHECK(a->equals(b));
        CHECK_FALSE(a->equals(c));
    }

    SECTION("add_to method") {
        auto a = NEW(NumVal)(10);
        auto b = NEW(NumVal)(5);
        auto result = std::dynamic_pointer_cast<NumVal>(a->add_to(b));
        CHECK(result->numVal == 15);
    }

    SECTION("mult_with method") {
        auto a = NEW(NumVal)(10);
        auto b = NEW(NumVal)(2);
        auto result = std::dynamic_pointer_cast<NumVal>(a->mult_with(b));
        CHECK(result->numVal == 20);
    }

    SECTION("print method") {
        auto a = NEW(NumVal)(10);
        std::stringstream ss;
        a->print(ss);
        CHECK(ss.str() == "10");
    }

    SECTION("is_true method") {
        auto a = NEW(NumVal)(10);
        CHECK_THROWS_AS(a->is_true(), std::runtime_error);
    }
}


TEST_CASE("BoolVal") {
    SECTION("constructor/print") {
        auto trueVal = NEW(BoolVal)(true);
        auto falseVal = NEW(BoolVal)(false);
        std::ostringstream outputTrue, outputFalse;

        trueVal->print(outputTrue);
        falseVal->print(outputFalse);

        CHECK(outputTrue.str() == "1");
        CHECK(outputFalse.str() == "0");
    }

    SECTION("equals") {
        CHECK(NEW(BoolVal)(true)->equals(NEW(BoolVal)(true)));
        CHECK_FALSE(NEW(BoolVal)(true)->equals(NEW(BoolVal)(false)));
    }

    SECTION("add_to") {
        auto boolVal = NEW(BoolVal)(true);
        auto anotherBoolVal = NEW(BoolVal)(false);

        CHECK_THROWS_AS(boolVal->add_to(anotherBoolVal), std::runtime_error);
    }

    SECTION("mult_with") {
        auto boolVal = NEW(BoolVal)(true);
        CHECK_THROWS_AS(boolVal->mult_with(NEW(BoolVal)(false)), std::runtime_error);
    }
}

TEST_CASE("Testing BoolExpr") {
    SECTION("constructor, print") {
        auto trueExpr = NEW(BoolExpr)(true);
        auto falseExpr = NEW(BoolExpr)(false);
        std::ostringstream outputTrue, outputFalse;
        trueExpr->print(outputTrue);
        falseExpr->print(outputFalse);
        CHECK(outputTrue.str() == "_true");
        CHECK(outputFalse.str() == "_false");
    }

    SECTION("equals") {
        CHECK(NEW(BoolExpr)(true)->equals(NEW(BoolExpr)(true)));
        CHECK_FALSE(NEW(BoolExpr)(true)->equals(NEW(BoolExpr)(false)));
    }

    SECTION("interp") {
        auto boolExpr = NEW(BoolExpr)(true);
        auto val = boolExpr->interp(Env::empty);
        auto boolVal = std::dynamic_pointer_cast<BoolVal>(val);
        REQUIRE(boolVal != nullptr);
        std::ostringstream output;
        boolVal->print(output);
        REQUIRE(output.str() == "1");
    }
}

TEST_CASE("Testing EqExpr") {
    auto varX = NEW(Var)("x");
    auto num1 = NEW(Num)(1);
    auto eqExpr = NEW(EqExpr)(varX, num1);

    SECTION("constructor and print") {
        std::ostringstream output;
        eqExpr->print(output);
        CHECK(output.str() == "x == 1");
    }

    SECTION("equals") {
        auto similarExpr = NEW(EqExpr)(varX, num1);
        auto differentExpr = NEW(EqExpr)(varX, NEW(Num)(2));
        CHECK(eqExpr->equals(similarExpr));
        CHECK_FALSE(eqExpr->equals(differentExpr));
    }

    SECTION("interp") {
        CHECK_THROWS_AS(eqExpr->interp(Env::empty), std::runtime_error);
    }
}

TEST_CASE("Testing IfExpr") {
    auto varX = NEW(Var)("x");
    auto num1 = NEW(Num)(1);
    auto num2 = NEW(Num)(2);
    auto condition = NEW(EqExpr)(varX, num1);
    auto ifExpr = NEW(IfExpr)(condition, num1, num2);

    SECTION("constructor, print") {
        std::ostringstream output;
        ifExpr->print(output);
        CHECK(output.str() == "_if x == 1 _then 1 _else 2");
    }

    SECTION("equals") {
        auto similarExpr = NEW(IfExpr)(condition, num1, num2);
        auto differentExpr = NEW(IfExpr)(NEW(EqExpr)(varX, num2), num1, num2);
        CHECK(ifExpr->equals(similarExpr));
        CHECK_FALSE(ifExpr->equals(differentExpr));
    }

    SECTION("interp") {
        CHECK_THROWS_AS(ifExpr->interp(Env::empty), std::runtime_error);
    }

}


TEST_CASE("Parse If, Bool and Functions") {

    SECTION("Parsing IfExpr") {
        CHECK(parse_str("_if _true _then 4 _else 5")->equals(
                NEW(IfExpr)(NEW(BoolExpr)(true), NEW(Num)(4), NEW(Num)(5))));
        CHECK(parse_str("_if _false _then 4 _else 5")->equals(
                NEW(IfExpr)(NEW(BoolExpr)(false), NEW(Num)(4), NEW(Num)(5))));
        CHECK(parse_str("_true")->equals(NEW(BoolExpr)(true)));
        CHECK(parse_str("_false")->equals(NEW(BoolExpr)(false)));
    }

    SECTION("Parsing EqExpr") {
        CHECK(parse_str("1 == 2")->interp(Env::empty)->equals(NEW(BoolVal)(false)));
        CHECK(parse_str("2 == 2")->interp(Env::empty)->equals(NEW(BoolVal)(true)));
        CHECK(parse_str("_if 1 == 2 _then 3 _else 4")->interp(Env::empty)->to_string() == "4");
        CHECK(parse_str("1 + 2 == 3 + 0")->interp(Env::empty)->equals(NEW(BoolVal)(true)));

        CHECK(parse_str("_if 1==1 _then 1 _else 2")->interp(Env::empty)->equals(NEW(NumVal)(1)));
    }

    SECTION("If Expr Interp") {
        CHECK(parse_str("_if 1==1 _then 1 _else 2")->interp(Env::empty)->equals(NEW(NumVal)(1)));
        CHECK(parse_str("_if 10==12 _then 7 _else 5")->interp(Env::empty)->equals(NEW(NumVal)(5)));
        CHECK(parse_str("_if 0==0 _then 14 _else 7")->interp(Env::empty)->equals(NEW(NumVal)(14)));
        CHECK(parse_str("_if -4==-5 _then 6 _else 8")->interp(Env::empty)->equals(NEW(NumVal)(8)));
    }
}



TEST_CASE("Function Print tests") {
    auto varX = NEW(Var)("x");
    auto num1 = NEW(Num)(1);
    auto body = NEW(Add)(varX, num1);
    auto body1 = NEW(Mult)(varX, num1);
    auto funExpr = NEW(FunExpr)("x", body);
    auto fun1Expr = NEW(FunExpr)("y", body1);

    SECTION("constructor, print") {
        std::ostringstream output;
        funExpr->print(output);
        CHECK(output.str() == "(_fun (x) (x+1))");

        output.str(""); // Clearing the output stream to test printing of the second function
        fun1Expr->print(output);
        CHECK(output.str() == "(_fun (y) (x*1))");
    }
}


TEST_CASE("Professor's given tests") {
    SECTION("Given Tests Assignment 2") {
        CHECK(NEW(Var)("x")->equals(NEW(Var)("x")));
        CHECK_FALSE(NEW(Var)("x")->equals(NEW(Var)("y")));
        CHECK_FALSE(NEW(Num)(1)->equals(NEW(Var)("x")));
        CHECK(NEW(Add)(NEW(Num)(2), NEW(Num)(3))->equals(NEW(Add)(NEW(Num)(2), NEW(Num)(3))));
        CHECK(NEW(Add)(NEW(Num)(2), NEW(Num)(3))->equals(NEW(Add)(NEW(Num)(3), NEW(Num)(2))));
    }

    SECTION("Given Tests Assignment 3") {
        CHECK(NEW(Mult)(NEW(Num)(3), NEW(Num)(2))->interp(Env::empty)->equals(NEW(NumVal)(6)));
        CHECK(NEW(Add)(NEW(Add)(NEW(Num)(10), NEW(Num)(15)), NEW(Add)(NEW(Num)(20), NEW(Num)(20)))->interp(Env::empty)->equals(NEW(NumVal)(65)));
        CHECK_THROWS_WITH(NEW(Var)("x")->interp(Env::empty), "no value for variable");
    }

    SECTION("Given Tests for Assignment 4") {
        CHECK(NEW(Num)(10)->to_string() == "10");
        std::stringstream ss1;
        NEW(Add)(NEW(Num)(1), NEW(Mult)(NEW(Num)(2), NEW(Num)(3)))->pretty_print_at(ss1);
        CHECK(ss1.str() == "1 + 2 * 3");

        std::stringstream ss2;
        NEW(Mult)(NEW(Num)(1), NEW(Add)(NEW(Num)(2), NEW(Num)(3)))->pretty_print_at(ss2);
        CHECK(ss2.str() == "1 * (2 + 3)");

        CHECK(NEW(Mult)(NEW(Num)(1), NEW(Add)(NEW(Num)(2), NEW(Num)(3)))->to_pp_string() == "1 * (2 + 3)");
        CHECK(NEW(Mult)(NEW(Mult)(NEW(Num)(8), NEW(Num)(1)), NEW(Var)("y"))->to_pp_string() == "(8 * 1) * y");
        CHECK(NEW(Mult)(NEW(Add)(NEW(Num)(3), NEW(Num)(5)), NEW(Mult)(NEW(Num)(6), NEW(Num)(1)))->to_pp_string() == "(3 + 5) * 6 * 1");
        CHECK(NEW(Mult)(NEW(Mult)(NEW(Num)(7), NEW(Num)(7)), NEW(Add)(NEW(Num)(9), NEW(Num)(2)))->to_pp_string() == "(7 * 7) * (9 + 2)");
    }

    SECTION("Nabil Given Test Assignment 5") {
        auto expr = NEW(Mult)(
                NEW(Mult)(NEW(Num)(2), NEW(_Let)("x", NEW(Num)(5), NEW(Add)(NEW(Var)("x"), NEW(Num)(1)))),
                NEW(Num)(3)
        );
        CHECK(expr->to_pp_string() == "(2 * _let x = 5\n      _in  x + 1) * 3");
    }

    SECTION("HW9 tests from Assignment Description") {
        std::string input = "1==2+3";
        auto expr = parse_str(input);
        CHECK(expr->interp(Env::empty)->to_string() == "0"); // False!

        std::string inputTwo = "1+1 == 2+0";
        auto exprTwo = parse_str(inputTwo);
        CHECK(exprTwo->interp(Env::empty)->to_string() == "1"); // True!
    }
}
