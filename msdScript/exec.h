//
// Created by Samantha Pope on 2/21/24.
//

#ifndef NEWMSDSCRIPT_EXEC_H
#define NEWMSDSCRIPT_EXEC_H

#include <string>
#include "pointer.h"

class ExecResult {
public:
    int exit_code;
    std::string out;
    std::string err;
    ExecResult() {
        exit_code = 0;
        out = "";
        err = "";
    }
};

extern ExecResult exec_program(int argc, const char * const *argv, std::string input);
ExecResult exec_program_wrapper(const std::string& path, const std::vector<std::string>& args, const std::string& input);

#endif //NEWMSDSCRIPT_EXEC_H
