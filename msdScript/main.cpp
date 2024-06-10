#define CATCH_CONFIG_RUNNER

#include <iostream>
#include <sstream> // Make sure to include this for std::istringstream
#include "cmdline.h"
#include "parse.h"
#include "Expr.h"
#include "catch.h"
#include "val.h"
#include "pointer.h"

int main(int argc, const char *argv[]) {
    try {
        run_mode_t mode = useArgs(argc, argv);

        switch (mode) {
            case do_test:
                if (Catch::Session().run(1, argv) != 0) {
                    // If Catch returns a non-zero value, then exit with failure
                    return 1;
                }
                break;

            case do_interp:
            case do_print:
            case do_pretty_print: {
                std::string line;
                if (std::getline(std::cin, line)) {
                    std::istringstream expr_stream(line);
                    PTR(Expr) expr = parse_expr(expr_stream);
                    if (mode == do_interp) {
                        PTR(Val) result = expr->interp(Env::empty); // Interpret the expression
                        std::string resultStr = result->to_string(); // Convert the result to a string
                        std::cout << resultStr << std::endl; // Print the result string
                    } else if (mode == do_print) {
                        expr->print(std::cout);
                        std::cout << std::endl;
                    } else { //pretty print method
                        std::cout << expr->to_pp_string() << std::endl;
                    }
                } else {
                    std::cerr << "No input received.\n";
                    return 1;
                }
                break;
            }

            default:
                std::cerr << "error: no valid mode selected \n";
                return 1;
        }
        return 0;
    } catch (const std::runtime_error& exn) {
        std::cerr << exn.what() << "\n";
        return 1;
    }
}
