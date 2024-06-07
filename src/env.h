#pragma once
#include <raisim/RaisimServer.hpp>
#include <raisim/World.hpp>

#include "controller.h"

namespace demo {
class DemoEnv {
 public:
  DemoEnv(std::string& urdf_path);
  virtual ~DemoEnv();
  void main();

 protected:
  Controller controller_;
  raisim::World* world_{nullptr};
  raisim::RaisimServer* server_{nullptr};
  ////////////////////////////////
  raisim::ArticulatedSystem* robot_{nullptr};
  Eigen::VectorXd gc_;
  Eigen::VectorXd gv_;
  Eigen::VectorXd gf_;
  int num_joints_{1};
};

}  // namespace demo
