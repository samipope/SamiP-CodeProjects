# Define the compiler to use
CXX = g++

# Compiler flags add
# -std=c++11 : Compile for C++11 standard
# -Wall : Turn on most compiler warnings
# -fsanitize=undefined : Check for undefined behavior
# -fno-sanitize-recover=undefined : Make the program exit immediately on undefined behavior
CXXFLAGS = -std=c++11 -Wall -fsanitize=undefined -fno-sanitize-recover=undefined

# Define any directories containing header files
INCLUDES =

# Define source (.cpp) and header (.hpp) files for the main target
SRCS_MAIN = main.cpp cmdline.cpp Expr.cpp ExprTest.cpp parse.cpp exec.cpp val.cpp Env.cpp # Add exec.cpp if it exists
HEADERS_MAIN = cmdline.h Expr.h catch.h parse.h exec.h val.h Env.h # Add exec.h if it exists

# Define source files for the test_msdscript target
SRCS_TEST_MSDSCRIPT = test_msdscript.cpp parse.cpp exec.cpp Expr.cpp val.cpp Env.cpp# Add any source files test_msdscript depends on
HEADERS_TEST_MSDSCRIPT = parse.h exec.h Expr.h val.h Env.h# Add any header files test_msdscript depends on

# Define the target executables
TARGET_MAIN = msdscript
TARGET_TEST_MSDSCRIPT = test_msdscript

# Define the part of the makefile that compiles the executable
all: $(TARGET_MAIN) $(TARGET_TEST_MSDSCRIPT)

$(TARGET_MAIN): $(SRCS_MAIN) $(HEADERS_MAIN)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET_MAIN) $(SRCS_MAIN)

$(TARGET_TEST_MSDSCRIPT): $(SRCS_TEST_MSDSCRIPT) $(HEADERS_TEST_MSDSCRIPT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET_TEST_MSDSCRIPT) $(SRCS_TEST_MSDSCRIPT)

test: $(TARGET_MAIN)
	./$(TARGET_MAIN) --test

# Target for generating documentation with Doxygen
doc:
	doxygen /Users/samanthapope/msdscriptRepo/msdScriptDocs/msdScriptDoxyFile

# Define clean command
.PHONY: clean
clean:
	rm -f $(TARGET_MAIN) $(TARGET_TEST_MSDSCRIPT) *.o

# Define a phony target to avoid problems with files named "all" or "clean"
.PHONY: all clean
