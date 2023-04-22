#include <iostream>
#include <fstream>
#include "ArgumentReader/ArgumentReader.h"

using namespace HwaUtil;
using namespace std;

int main() {
    ArgumentReader argReader;
    argReader.AddArg("calculation");
    argReader.AddArg("matrix_type");
    argReader.AddArg("matrix_print");
    argReader.AddArg("nrows");
    argReader.AddArg("ncols");


    ifstream fs("matmul5.txt");

    argReader.ReadArgs(fs);

    fs.close();

    std::cout << argReader << std::endl;

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
