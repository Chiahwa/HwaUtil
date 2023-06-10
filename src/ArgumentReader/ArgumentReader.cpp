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

    int ArgumentReader::ReadArgs(istream &is) {
        Timer::tick("HwaUtil::ArgumentReader", "ReadArgs");
        constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

        if (NArgs == 0) {
            Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
            return 0;
        }

        bool data_label_found = !DataLabel.empty();

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

            if(data_label_found && name.rfind(DataLabel, 0) == 0)
                break;

            if (!ArgID.contains(name)) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::out_of_range("Argument " + name + " does not exist!");
            }
            if (ArgVal.contains(ArgID[name])) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Argument " + name + " declared more than once!");
            }

            string val;
            ssline >> val;
            if (ssline.fail()) {
                Timer::tock("HwaUtil::ArgumentReader", "ReadArgs");
                throw std::runtime_error("Value of argument " + name + " missing!");
            }
            /* data beginning with '$' is case-sensitive */
            if (val[0] == '$')
                val.erase(val.begin());
            else
                transform(val.begin(), val.end(), val.begin(), ::tolower);

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
        return nline;
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

    ArgumentReader::ArgumentReader() = default;

    bool ArgumentReader::SetDataLabel(string label) {
        Timer::tick("HwaUtil::ArgumentReader", "SetDataLabel");
        transform(label.begin(),label.end(), label.begin(), ::tolower);
        if (!label.empty()) {
            label.erase(0, label.find_first_not_of(' '));
            label.erase(label.find_last_not_of(' ') + 1);
        }
        if(ArgID.contains(label)) {
            Timer::tock("HwaUtil::ArgumentReader", "SetDataLabel");
            return false;
        }else {
            DataLabel=label;
            Timer::tock("HwaUtil::ArgumentReader", "SetDataLabel");
            return true;
        }
    }


} // HwaUtil