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
        if (strcmp(argv[i], "--test") == 0) {
            mode = do_test;
        }

        else if(strcmp(argv[i], "--help") == 0){
            mode = do_help;
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
