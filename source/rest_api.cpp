#include "rest_api.h"
#include <iostream>
#include <string>
#include <sstream>
#include <string.h>

extern CfgParameter CfgParmsDef[];
extern void process_drc(); // We'll need to define this in drc.cpp

std::string get_config_json() {
    std::stringstream ss;
    ss << "{";
    int i = 0;
    bool first = true;
    while (CfgParmsDef[i].PType != CfgEnd) {
        if (!first) ss << ", ";
        ss << "\"" << CfgParmsDef[i].PName << "\": ";
        switch (CfgParmsDef[i].PType) {
            case CfgInt:
                ss << (*((int*)CfgParmsDef[i].PValue));
                break;
            case CfgFloat:
                ss << (*((float*)CfgParmsDef[i].PValue));
                break;
            case CfgDouble:
                ss << (*((double*)CfgParmsDef[i].PValue));
                break;
            case CfgString: {
                char* str = *((char**)CfgParmsDef[i].PValue);
                ss << "\"" << (str ? str : "") << "\"";
                break;
            }
            default:
                ss << "null";
                break;
        }
        first = false;
        i++;
    }
    ss << "}";
    return ss.str();
}

bool update_config(const std::string& key, const std::string& value) {
    int i = 0;
    while (CfgParmsDef[i].PType != CfgEnd) {
        if (std::string(CfgParmsDef[i].PName) == key) {
            switch (CfgParmsDef[i].PType) {
                case CfgInt:
                    *((int*)CfgParmsDef[i].PValue) = std::stoi(value);
                    return true;
                case CfgFloat:
                    *((float*)CfgParmsDef[i].PValue) = std::stof(value);
                    return true;
                case CfgDouble:
                    *((double*)CfgParmsDef[i].PValue) = std::stod(value);
                    return true;
                case CfgString: {
                    char** str_ptr = (char**)CfgParmsDef[i].PValue;
                    if (*str_ptr != NULL) free(*str_ptr);
                    *str_ptr = strdup(value.c_str());
                    return true;
                }
                default:
                    return false;
            }
        }
        i++;
    }
    return false;
}

void start_rest_server(int port) {
    httplib::Server svr;

    svr.Get("/config", [](const httplib::Request& req, httplib::Response& res) {
        res.set_content(get_config_json(), "application/json");
    });

    // Handle updates like POST /config?MPLowerWindow=123
    svr.Post("/config", [](const httplib::Request& req, httplib::Response& res) {
        for (auto it = req.params.begin(); it != req.params.end(); ++it) {
            if (update_config(it->first, it->second)) {
                std::cout << "Updated " << it->first << " to " << it->second << std::endl;
            } else {
                std::cout << "Failed to update " << it->first << std::endl;
            }
        }
        res.set_content(get_config_json(), "application/json");
    });

    svr.Post("/update_output", [](const httplib::Request& req, httplib::Response& res) {
        std::cout << "Triggering recalculation..." << std::endl;
        process_drc();
        res.set_content("{\"status\": \"success\"}", "application/json");
    });

    std::cout << "Starting REST server on port " << port << "..." << std::endl;
    svr.listen("0.0.0.0", port);
}
