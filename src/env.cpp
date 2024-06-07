#include "env.h"

#include "pinocchio/algorithm/joint-configuration.hpp"

namespace demo {

DemoEnv::DemoEnv(std::string& urdf_path) {
  controller_.init(urdf_path);
  world_ = new raisim::World();
  ////////////////////////////////
  robot_ = world_->addArticulatedSystem(urdf_path);
  robot_->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
  robot_->setGeneralizedForce(Eigen::VectorXd::Zero(robot_->getDOF()));
  num_joints_ = robot_->getDOF();
  gc_.setZero(num_joints_);
  gv_.setZero(num_joints_);
  gf_.setZero(num_joints_);
  robot_->getState(gc_, gv_);
  ///////////////////////////////
  world_->addGround();
  world_->setTimeStep(0.005);
  raisim::Vec<3> gravity{0, 0, -9.81};
  world_->setGravity(gravity);
  ///////////////////////////////
  server_ = new raisim::RaisimServer(world_);
}

DemoEnv::~DemoEnv() {
  if (server_) {
    server_->killServer();
    delete server_;
  }
  delete world_;
}

void DemoEnv::main() {
  Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
  Eigen::Vector3d translation{0.3, 0., 0.9};
  raisim::Visuals* target = server_->addVisualSphere("target", 0.05);
  target->setPosition(translation);
  target->setColor(1.0, 0., 0., 0.92);
  pinocchio::SE3 desired_SE3(rot, translation);
  Eigen::VectorXd origin_q = pinocchio::neutral(controller_.getRobotModel());
  Eigen::VectorXd desired_theta =
      controller_.solveInvKinematic(desired_SE3, origin_q);
  server_->launchServer();
  raisim::MSLEEP(1000);
  while (1) {
    // robot_->setGeneralizedCoordinate(desired_theta);
    raisim::MSLEEP(4);
    gf_.tail(num_joints_) =
        controller_.calculateTorque(desired_theta, gc_, gv_);
    robot_->setGeneralizedForce(gf_);
    if (server_) server_->lockVisualizationServerMutex();
    world_->integrate();
    if (server_) server_->unlockVisualizationServerMutex();
    robot_->getState(gc_, gv_);
  }
}

}  // namespace demo
