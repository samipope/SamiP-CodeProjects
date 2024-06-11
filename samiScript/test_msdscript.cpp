#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <sstream>
#include "exec.h"

std::string random_expr_string(int depth = 0);
std::string random_variable();

std::string random_expr_string(int depth) {
    int r = rand() % 10; // Adjust ratios as needed
    if (depth > 2 || r < 3) { // Base case: favor numbers at deeper recursion levels
        return std::to_string(rand() % 100); // Numbers from 0 to 99
    } else if (r < 6) { // Recursive case: addition
        return random_expr_string(depth + 1) + " + " + random_expr_string(depth + 1);
    } else if (r < 8) { // Variable
        return random_variable();
    } else { // Let expression
        std::string var = random_variable();
        std::string valExpr = random_expr_string(depth + 1);
        std::string bodyExpr = random_expr_string(depth + 1);
        return "let " + var + " = " + valExpr + " in " + bodyExpr;
    }
}

std::string random_variable() {
    // Generates a single lowercase letter as a variable
    return std::string(1, 'a' + (rand() % 26));
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <path_to_msdscript1> <path_to_msdscript2>" << std::endl;
        return 1;
    }

    srand(static_cast<unsigned int>(time(nullptr))); // Initialize random seed

    // Use the paths provided by the user
    std::string msdscriptPath = argv[1]; // Path to the first msdscript executable
    std::string msdscript2Path = argv[2];
    for (int i = 0; i < 100; ++i) {
        std::string expr = random_expr_string();
        std::cout << "Testing expression: " << expr << "\n";

        ExecResult result2 = exec_program_wrapper(msdscript2Path, {"--interp"}, expr + "\n");
        ExecResult result1 = exec_program_wrapper(msdscriptPath, {"--interp"}, expr + "\n");


        if (result1.out != result2.out) {
            std::cerr << "Mismatch found!\nExpression: " << expr << "\nResult1: " << result1.out << "\nResult2: "
                      << result2.out << std::endl;}
//        } else {
//            std::cout << "Match: " << expr << " => " << result1.out << std::endl;
//        }

    }

    return 0;
}
