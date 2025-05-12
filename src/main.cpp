#include "amflow.hpp"


int main(int argc, const char** argv) {
    if (argc != 2)
        amflow::usage();

    amflow solver(argv[1]);
    
}

