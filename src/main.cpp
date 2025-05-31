#include "amflow.hpp"


int main(int argc, const char** argv) {
    if (argc != 2)
        amflow::usage();

    amflow solver(argv[1], {{1,1,1,1,1,1,1,-3,0},{1,1,1,1,1,1,1,-2,-1},{1,1,1,1,1,1,1,-1,-2},{1,1,1,1,1,1,1,0,-3}}, {}, "/root/amflow_cpp_test5");
    solver.forward();
    solver.backward(GiNaC::numeric(1)/1000);
}

