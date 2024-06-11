#define CATCH_CONFIG_RUNNER

#include <iostream>
#include <sstream> // Make sure to include this for std::istringstream
#include "cmdline.h"
#include "parse.h"
#include "Expr.h"
#include "catch.h"
#include "val.h"
#include "pointer.h"

using namespace std;
//
//int main(int argc, const char *argv[]) {
//    try {
//        run_mode_t mode = useArgs(argc, argv);
//
//        switch (mode) {
//            case do_help: {
//                cout << "--test: tests the code.\n";
//                cout << "--help: check your options.\n";
//                cout << "--print: prints your input out.\n";
//                cout << "--pretty_print: prints in a prettier way, helpful for large functions.\n";
//                break;
//            }
//
//            case do_test: {
//                if (Catch::Session().run(1, argv) != 0) {
//                    // If Catch returns a non-zero value, then exit with failure
//                    return 1;
//                }
//                break;
//            }
//            case do_interp: {
//                std::string line;
//                if (std::getline(std::cin, line)) {
//                    std::istringstream expr_stream(line);
//                    PTR(Expr) expr = parse_expr(expr_stream);
//                PTR(Val) result = expr->interp(Env::empty); // Interpret the expression
//                std::string resultStr = result->to_string(); // Convert the result to a string
//                std::cout << resultStr << std::endl; // Print the result string
//                PTR(Expr) e = parse(std::cin);
//                cout << e->interp(Env::empty)->to_string() << "\n";
//                break;
//            }
//            case do_print: {
//                PTR(Expr) e = parse(std::cin);
//                std::cout << e->to_string() << "\n";
//                break;
//            }
//            case do_pretty_print: {
//                std::string line;
//                if (std::getline(std::cin, line)) {
//                    std::istringstream expr_stream(line);
//                    PTR(Expr) expr = parse_expr(expr_stream);
//                    if (mode == do_interp) {
//                        PTR(Val) result = expr->interp(Env::empty); // Interpret the expression
//                        std::string resultStr = result->to_string(); // Convert the result to a string
//                        std::cout << resultStr << std::endl; // Print the result string
//                    } else if (mode == do_print) {
//                        expr->print(std::cout);
//                        std::cout << std::endl;
//                    } else { //pretty print method
//                        std::cout << expr->to_pp_string() << std::endl;
//                    }
//                } else {
//                    std::cerr << "No input received.\n";
//                    return 1;
//                }
//                break;
//            }
//
//            default:
//                std::cerr << "error: no valid mode selected \n";
//                return 1;
//        }
//        return 0;
//    } catch (const std::runtime_error& exn) {
//        std::cerr << exn.what() << "\n";
//        return 1;
//    }
//}

    void handleHelp() {
        cout << "--test: tests the code.\n";
        cout << "--help: check your options.\n";
        cout << "--print: prints your input out.\n";
        cout << "--pretty_print: prints in a prettier way, helpful for large functions.\n";
    }

    int handleTest(const char *argv[]) {
        return Catch::Session().run(1, argv) != 0 ? 1 : 0;
    }

    int handleInterp() {
        std::string line;
        if (std::getline(std::cin, line)) {
            std::istringstream expr_stream(line);
            PTR(Expr) expr = parse_expr(expr_stream);
//            std::cout << "here is the xpr that was parsed" + expr->to_string();
            PTR(Val) result = expr->interp(Env::empty);
//            std::cout << "here is the result: ";
            std::cout << result->to_string() << std::endl;
        } else {
            std::cerr << "No input received.\n";
            return 1;
        }
        return 0;
    }

    int handlePrint() {
        PTR(Expr) e = parse(std::cin);
        std::cout << e->to_string() << "\n";
        return 0;
    }

    int handlePrettyPrint() {
        std::string line;
        if (std::getline(std::cin, line)) {
            std::istringstream expr_stream(line);
            PTR(Expr) expr = parse_expr(expr_stream);
            std::cout << expr->to_pp_string() << std::endl;
        } else {
            std::cerr << "No input received.\n";
            return 1;
        }
        return 0;
    }

    int main(int argc, const char *argv[]) {
        try {
            run_mode_t mode = useArgs(argc, argv);

            switch (mode) {
                case do_help:
                    handleHelp();
                    break;
                case do_test:
                    return handleTest(argv);
                case do_interp:
                    return handleInterp();
                case do_print:
                    return handlePrint();
                case do_pretty_print:
                    return handlePrettyPrint();
                default:
                    cerr << "error: no valid mode selected \n";
                    return 1;
            }
            return 0;
        } catch (const std::runtime_error& exn) {
            cerr << exn.what() << "\n";
            return 1;
        }
    }