#include "amflow.hpp"


int main(int argc, const char** argv) {
    if (argc != 2)
        amflow::usage();

    amflow solver(argv[1], {{1,1,1,2,1,1,1},{1,2,1,1,1,1,0}}, {}, "/root/amflow_cpp_test3");
    solver.forward();
}

