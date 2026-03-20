#include "rest_api.h"
#include <iostream>
#include <string>
#include <sstream>
#include <string.h>

extern CfgParameter CfgParmsDef[];
extern void trigger_recalculation(); // We'll need to define this in drc.cpp

// Important parameters (marked with * in docs) and their valid min/max ranges based on documentation
struct ConfigMeta {
    const char* name;
    bool significant;
    float min_val;
    float max_val;
};

ConfigMeta config_metadata[] = {
    {"BCImpulseCenter", true, 0, -1}, // Variable, depends on file
    {"MPLowerWindow", true, 16384, 65536},
    {"MPUpperWindow", true, 22, 128},
    {"MPWindowExponent", true, 0.7, 1.2},
    {"MPFSharpness", true, 0.1, 0.75},
    {"EPLowerWindow", true, 1024, 4096},
    {"EPUpperWindow", true, 22, 128},
    {"EPWindowExponent", true, 0.5, 1.2},
    {"EPFSharpness", true, 0.1, 0.75},
    {"PTReferenceWindow", true, 100, 65536}, // Doc says default 26460, typical 150-500ms
    {"PTBandWidth", true, -3, 1}, // -2 is ERB, -1 is Bark, >0 fractions of octaves
    {"PTPeakDetectionStrength", true, 5, 30}, // "Values above 50 are probably going to cause numerical problems"
    {"RTLowerWindow", true, 16384, 65536}, // Same as MPLowerWindow
    {"RTUpperWindow", true, 22, 128}, // Same as MPUpperWindow
    {"RTWindowExponent", true, 0.7, 1.2}, // Same as MPWindowExponent
    {"RTFSharpness", true, 0.1, 0.75}, // Same as MPFSharpness
    {"PSPointsFile", true, 0, 0}, // File path
    {NULL, false, 0, 0}
};

const ConfigMeta* get_meta(const std::string& name) {
    for (int i = 0; config_metadata[i].name != NULL; ++i) {
        if (name == config_metadata[i].name) {
            return &config_metadata[i];
        }
    }
    return NULL;
}

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

std::string get_schema_json() {
    std::stringstream ss;
    ss << "[\n";
    int i = 0;
    bool first = true;
    while (CfgParmsDef[i].PType != CfgEnd) {
        if (!first) ss << ",\n";
        
        std::string name = CfgParmsDef[i].PName;
        const ConfigMeta* meta = get_meta(name);
        
        ss << "  {\n";
        ss << "    \"name\": \"" << name << "\",\n";
        
        // Type
        ss << "    \"type\": ";
        switch (CfgParmsDef[i].PType) {
            case CfgInt: ss << "\"int\",\n"; break;
            case CfgFloat: ss << "\"float\",\n"; break;
            case CfgDouble: ss << "\"double\",\n"; break;
            case CfgString: ss << "\"string\",\n"; break;
            default: ss << "\"unknown\",\n"; break;
        }
        
        // Meta info
        if (meta && meta->significant) {
            ss << "    \"significant\": true";
            if (meta->min_val != meta->max_val) {
                ss << ",\n    \"min\": " << meta->min_val << ",\n";
                ss << "    \"max\": " << meta->max_val << "\n";
            } else {
                ss << ",\n    \"min\": null,\n";
                ss << "    \"max\": null\n";
            }
        } else {
            ss << "    \"significant\": false,\n";
            ss << "    \"min\": null,\n";
            ss << "    \"max\": null\n";
        }
        
        ss << "  }";
        first = false;
        i++;
    }
    ss << "\n]";
    return ss.str();
}

bool update_config(const std::string& key, const std::string& value) {
    int i = 0;
    while (CfgParmsDef[i].PType != CfgEnd) {
        if (std::string(CfgParmsDef[i].PName) == key) {
            const ConfigMeta* meta = get_meta(key);
            switch (CfgParmsDef[i].PType) {
                case CfgInt: {
                    int val = std::stoi(value);
                    if (meta && meta->significant && meta->min_val != meta->max_val) {
                        if (val < meta->min_val) val = meta->min_val;
                        if (val > meta->max_val) val = meta->max_val;
                    }
                    *((int*)CfgParmsDef[i].PValue) = val;
                    return true;
                }
                case CfgFloat: {
                    float val = std::stof(value);
                    if (meta && meta->significant && meta->min_val != meta->max_val) {
                        if (val < meta->min_val) val = meta->min_val;
                        if (val > meta->max_val) val = meta->max_val;
                    }
                    *((float*)CfgParmsDef[i].PValue) = val;
                    return true;
                }
                case CfgDouble: {
                    double val = std::stod(value);
                    if (meta && meta->significant && meta->min_val != meta->max_val) {
                        if (val < meta->min_val) val = meta->min_val;
                        if (val > meta->max_val) val = meta->max_val;
                    }
                    *((double*)CfgParmsDef[i].PValue) = val;
                    return true;
                }
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
    
    svr.Get("/schema", [](const httplib::Request& req, httplib::Response& res) {
        res.set_content(get_schema_json(), "application/json");
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
        trigger_recalculation();
        res.set_content("{\"status\": \"success\"}", "application/json");
    });

    std::cout << "Starting REST server on port " << port << "..." << std::endl;
    svr.listen("0.0.0.0", port);
}
