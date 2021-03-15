//
// Created by 徐溶延 on 2020/11/17.
//

#include "parse_config_file.h"

#include <fstream>
#include <iostream>

#include "../common/string_utils.h"

#include <boost/program_options.hpp>

using namespace boost;

int parseConfigFile(const std::string &filename, common::argument &args) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cout << "unable to read the configure file" << std::endl;
        return __LINE__;
    }
    program_options::options_description config("command parameters");
    config.add_options()
            ("help, h", "hlep message")
            ("in_path, s", program_options::value<std::string>(), "input file path.")
            ("out_path, s", program_options::value<std::string>(), "output polyline file path.")
            ("data_path, s", program_options::value<std::string>(), "control points file")
            ("topo_path, s", program_options::value<std::string>(), "script of controlling topological operations")
            ("loop, i", program_options::value<int>()->default_value(1), ".");

    program_options::variables_map vm;
    program_options::store(program_options::parse_config_file(ifs, config), vm);
    ifs.close();
    program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << config << std::endl;
        return __LINE__;
    }
    args.in_path = DATA_SHARED_PATH"/" + vm["in_path"].as<std::string>();
    args.out_path = DATA_SHARED_PATH"/" + vm["out_path"].as<std::string>();
    args.data_path = DATA_SHARED_PATH"/" + vm["data_path"].as<std::string>();
//    args.topo_path = DATA_SHARED_PATH"/" + vm["topo_path"].as<std::string>();
//    args.loop = vm["loop"].as<int>();
    args.out_init_path = common::findPrefix(args.out_path) + "_init.vtk";
    return 0;
}
