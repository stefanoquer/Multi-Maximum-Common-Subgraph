#pragma once
#include <vector>
#include <string>
#include "mcsp.h"
#include "CLI11.hpp"

using namespace std;

class App : public CLI::App {
public:
    App(const std::string& desc, const std::string& filename);
    ~App() override;
    void run();
    args arguments;

private:
    bool _set_minmax = false;
    bool _set_minproduct = false;
    bool _set_minmin = false;
    bool _set_minsum = false;
};