//
// Created by Samantha Pope on 1/19/24.
//
#include "pointer.h"

#ifndef CS6015PROJECT_CMDLINE_H
#define CS6015PROJECT_CMDLINE_H

typedef enum {
    do_nothing,
    do_test,
    do_interp,
    do_print,
    do_pretty_print
} run_mode_t;

run_mode_t useArgs(int argc, const char * argv[]);


#endif //CS6015PROJECT_CMDLINE_H
