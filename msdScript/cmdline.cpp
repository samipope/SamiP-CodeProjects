//
// Created by Samantha Pope on 1/19/24.
//

#include <iostream>
#include <cstring> // For strcmp
#include "cmdline.h"
#include "catch.h"
#include "parse.h"
#include "pointer.h"

run_mode_t useArgs(int argc, const char *argv[]) {
    run_mode_t mode = do_nothing;
 //   static bool testAlreadyPressed = false;

    for (int i = 1; i < argc; ++i) {
        //loop for all of it
        if (strcmp(argv[i], "--help") == 0) {
            std::cout << "Usage: ./msdscript --[option]\n"
                      << "Options:\n"
                      << "  --help  Show this help message\n"
                      << "  --test  Run tests\n";
            exit(0);
        }

            //have to control for if it is the first test that they wrote or not
        else if (strcmp(argv[i], "--test") == 0) {
            mode = do_test;
        }

        else if(strcmp(argv[i], "--interp")==0){
            mode = do_interp;
        }

        else if(strcmp(argv[i], "--print")==0){
            mode = do_print;
        }

        else if (strcmp(argv[i],"--pretty_print")==0){
            mode = do_pretty_print;
        }

        else {
            std::cerr << "Error: Unknown argument:  " << argv[i] << "\n";
            exit(1);
        }

    }
    return mode;
}
