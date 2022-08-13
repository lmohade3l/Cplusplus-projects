#include <iostream>
#include<fstream>
#include "aphw1.h"
#include "gtest/gtest.h"


using Matrix = std::vector<std::vector<double>>;
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    std::cout << "RUNNING TESTS ..." << std::endl;
    int ret{RUN_ALL_TESTS()};
    if (!ret)
        std::cout << "<<<SUCCESS>>>" << std::endl;
    else
        std::cout << "FAILED" << std::endl;
    
    
    return 0;
}

