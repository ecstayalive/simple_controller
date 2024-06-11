#pragma once
#include <memory>
#include <raisim/RaisimServer.hpp>
#include <raisim/World.hpp>

#include "controller.h"
#include "kinematic.h"

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
  //////////////////////////////////////////////////////////////////////////
  raisim::ArticulatedSystem* robot_{nullptr};
  Eigen::VectorXd gc_;
  Eigen::VectorXd gv_;
  Eigen::VectorXd gf_;
  std::vector<std::string> arm_links_names_{"link01", "line02", "link03",
                                            "link04", "link05", "link06"};
  std::vector<size_t> arm_joints_frame_idx_;
  std::vector<size_t> arm_bodies_idx_;
  std::vector<Eigen::Vector3d> arm_joint_init_pos_;
  std::vector<Eigen::Matrix3d> arm_joint_init_rot_;
  std::vector<Eigen::Vector2d> arm_joint_limits_;
  Eigen::VectorXd arm_joint_init_angles_;
  std::vector<Eigen::Vector3d> arm_joint_axis_;
  Eigen::Vector3d end_effector_init_pos_P_;
  size_t end_effector_frame_idx;
  int num_joints_{1};
  //////////////////////////////////////////////////////////////////////////
  std::unique_ptr<KinematicChain> kinematic_chain_{nullptr};
};

}  // namespace demo
