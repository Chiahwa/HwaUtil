//
// Created by Chiahwa Young on 2023/3/15.
//

#include "ArgumentReader.h"
#include <exception>
#include <iostream>
#include <sstream>
#include <vector>


namespace HwaUtil {
    ArgumentReader::ArgumentReader() = default;

    //TODO：多余的参数检测
    void ArgumentReader::ReadArgs(istream &is) {
        constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

        if (NArgs == 0)return;
        string line;
        int nline=0;
        while (getline(is, line)) {
            ++nline;
            stringstream buf(line);
            getline(buf,line,'#');
            stringstream ssline(line);

            string name;
            ssline >> name;
            transform(name.begin(), name.end(), name.begin(), ::tolower);
            if (ssline.fail()) {
                //is.ignore(max_size, '\n');
                continue;
            }

            string val;
            ssline >> val;
            if (ssline.fail())
                throw std::runtime_error("Value of argument " + name + " missing!");
            transform(val.begin(), val.end(), val.begin(), ::tolower);

            if (!ArgID.contains(name))
                throw std::out_of_range("Argument " + name + " does not exist!");
            if (ArgVal.contains(ArgID[name]))
                throw std::runtime_error("Argument " + name + " declared more than once!");

            ArgVal[ArgID[name]] = val;
            //is.ignore(max_size, '\n');
            string other;
            ssline >> other;
            if(!other.empty()){
                throw std::runtime_error("Line " + to_string(nline) + ": too much input!");
            }
        }
        for (auto &arg: ArgID)
            if (!ArgVal.contains(arg.second))
                throw std::runtime_error("Argument " + arg.first + " missing!");
    }

    /* returns true if success */
    bool ArgumentReader::AddArg(const string &name) {
        if (ArgID.contains(name)) return false;
        else {
            ArgID[name] = ++NArgs;
            return true;
        }
    }

    string ArgumentReader::GetArgV(const string &name) {
        if (!ArgID.contains(name))
            throw std::out_of_range("GetArgV(name): Argument name does not exist!");
        return ArgVal[ArgID[name]];
    }

    string ArgumentReader::GetArgV(const int ID) {
        if (ID <= NArgs)
            return ArgVal[ID];
        else throw std::out_of_range("GetArgV(ID): Argument ID out of range!");
    }

} // HwaUtil