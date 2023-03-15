#include <iostream>
#include <fstream>
#include "src/ArgumentReader/ArgumentReader.h"
using namespace HwaUtil;
using namespace std;

int main() {
    ArgumentReader argReader;
    argReader.AddArg("calculation");
    argReader.AddArg("matrix_type");
    argReader.AddArg("matrix_print");
    argReader.AddArg("nrows");
    argReader.AddArg("ncols");

    ifstream fs("input.txt");

    argReader.ReadArgs(fs);

    fs.close();

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
