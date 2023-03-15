//
// Created by Chiahwa Young on 2023/3/15.
//

#ifndef HWAUTIL_ARGUMENTREADER_H
#define HWAUTIL_ARGUMENTREADER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
using namespace std;

namespace HwaUtil {

    class ArgumentReader {
    private:
        int NArgs=0;
        /* ArgID starts from 1 */
        map<string,int> ArgID;
        map<int,string> ArgVal;

    public:
        explicit ArgumentReader();
        bool AddArg(const string& name);
        void ReadArgs(istream &is);
        string GetArgV(const string &name);
        string GetArgV(const int ID);

    };

} // HwaUtil

#endif //HWAUTIL_ARGUMENTREADER_H
