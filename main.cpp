#include <iostream>

#include "src/env.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Require 1 argument: robot urdf file. Get " << argc - 1
              << ".\n";
    return -1;
  } else {
    std::cout << "Loading " << argv[1] << "!\n";
  }
  std::string urdf_path = argv[1];
  ////////////////////////////////////////////////////////////////////////////////////////////////
  demo::DemoEnv demo(urdf_path);
  demo.main();
  return 0;
}
