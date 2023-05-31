//
// Created by Chiahwa Young on 2023/3/15.
//

#include "ArgumentReader.h"
#include <exception>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Timer/Timer.h"


namespace HwaUtil {
    ArgumentReader::ArgumentReader() = default;


    void ArgumentReader::ReadArgs(istream &is) {
        Timer::tick("HwaUtil::ArgumentReader", "ReadArgs");
        constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

        if (NArgs == 0) {
            Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
            return;
        }
        string line;
        int nline = 0;
        while (getline(is, line)) {
            ++nline;
            stringstream buf(line);
            getline(buf, line, '#');
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
            if (ssline.fail()) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Value of argument " + name + " missing!");
            }
            if (val[0] == '$')
                val.erase(val.begin());
            else
                transform(val.begin(), val.end(), val.begin(), ::tolower);

            if (!ArgID.contains(name)) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::out_of_range("Argument " + name + " does not exist!");
            }
            if (ArgVal.contains(ArgID[name])) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Argument " + name + " declared more than once!");
            }

            ArgVal[ArgID[name]] = val;
            //is.ignore(max_size, '\n');
            string other;
            ssline >> other;
            if (!other.empty()) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Line " + to_string(nline) + ": too much input!");
            }
        }
        for (auto &arg: ArgID)
            if (!ArgVal.contains(arg.second)) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Argument " + arg.first + " missing!");
            }
        Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
    }

    /* returns true if success */
    bool ArgumentReader::AddArg(string name) {
        Timer::tick("HwaUtil::ArgumentReader", "AddArg");
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (ArgID.contains(name)) {
            Timer::tock("HwaUtil::ArgumentReader", "AddArg");
            return false;
        } else {
            ArgID[name] = ++NArgs;
            Timer::tock("HwaUtil::ArgumentReader", "AddArg");
            return true;
        }
    }

    string ArgumentReader::GetArgV(const string &name) {
        Timer::tick("HwaUtil::ArgumentReader", "GetArgV");
        if (!ArgID.contains(name)) {
            Timer::tock("HwaUtil::ArgumentReader", "GetArgV");
            throw std::out_of_range("GetArgV(name): Argument name does not exist!");
        }
        Timer::tock("HwaUtil::ArgumentReader", "GetArgV");
        return ArgVal[ArgID[name]];
    }

    string ArgumentReader::GetArgV(int ID) {
        Timer::tick("HwaUtil::ArgumentReader", "GetArgV");
        if (ID <= NArgs) {
            Timer::tock("HwaUtil::ArgumentReader", "GetArgV");
            return ArgVal[ID];
        } else {
            Timer::tock("HwaUtil::ArgumentReader", "GetArgV");
            throw std::out_of_range("GetArgV(ID): Argument ID out of range!");
        }
    }

    ostream &operator<<(ostream &os, const ArgumentReader &ar) {
        Timer::tick("HwaUtil::ArgumentReader", "operator<<");
        for (auto &arg: ar.ArgID) {
            os << arg.first << " = " << ar.ArgVal.at(arg.second) << endl;
        }
        Timer::tock("HwaUtil::ArgumentReader", "operator<<");
        return os;
    }


} // HwaUtil