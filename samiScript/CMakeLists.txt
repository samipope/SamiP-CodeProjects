cmake_minimum_required(VERSION 3.27)
project(newMSDSCRIPT)

set(CMAKE_CXX_STANDARD 17)

add_executable(newMSDSCRIPT main.cpp
        cmdline.cpp
        Expr.cpp
        parse.cpp
        ExprTest.cpp
        exec.h
        exec.cpp
        test_msdscript.cpp
        val.h
        val.cpp
        pointer.h
        Env.cpp
)
