#include "env.h"

#include <chrono>

#include "pinocchio/algorithm/frames.hpp"
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
  arm_joint_init_angles_.setOnes(num_joints_);
  ///////////////////////////////////////////////////////////////////////////////
  arm_joints_frame_idx_.reserve(6);
  for (auto link_name : arm_links_names_) {
    arm_joints_frame_idx_.push_back(robot_->getFrameIdxByName(link_name));
    arm_bodies_idx_.push_back(robot_->getBodyIdx(link_name));
  }
  end_effector_frame_idx = robot_->getFrameIdxByName("endEffector");
  std::vector<raisim::Vec<3>> rmt_joints_axis, rmt_joints_init_pos;
  std::vector<raisim::Mat<3, 3>> rmt_joints_init_rot;
  std::vector<raisim::Vec<2>> rmt_joints_limits;
  rmt_joints_axis = robot_->getJointAxis_P();
  rmt_joints_init_pos = robot_->getJointPos_P();
  rmt_joints_init_rot = robot_->getJointOrientation_P();
  for (int i = 1; i < rmt_joints_axis.size(); ++i) {
    arm_joint_axis_.emplace_back(rmt_joints_axis[i].e());
    arm_joint_init_pos_.emplace_back(rmt_joints_init_pos[i].e());
    arm_joint_init_rot_.emplace_back(rmt_joints_init_rot[i].e());
  }
  rmt_joints_limits = robot_->getJointLimits();
  for (int i = 0; i < rmt_joints_limits.size(); ++i)
    arm_joint_limits_.emplace_back(rmt_joints_limits[i].e());
  raisim::Vec<3> rmt_ee_pos_W, rmt_ee_pos_P;
  robot_->getFramePosition(end_effector_frame_idx, rmt_ee_pos_W);
  robot_->getPositionInBodyCoordinate(arm_bodies_idx_.back(), rmt_ee_pos_W,
                                      rmt_ee_pos_P);
  end_effector_init_pos_P_ = rmt_ee_pos_P.e();
  kinematic_chain_ = std::make_unique<KinematicChain>(
      arm_joint_axis_, arm_joint_init_rot_, arm_joint_init_pos_,
      end_effector_init_pos_P_, arm_joint_init_angles_, arm_joint_limits_);
  //////////////////////////////////////////////////////////////////////////////
  world_->addGround();
  world_->setTimeStep(0.005);
  raisim::Vec<3> gravity{0, 0, -9.81};
  world_->setGravity(gravity);
  //////////////////////////////////////////////////////////////////////////////
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
  Eigen::Quaterniond q = Eigen::Quaterniond(1, 0.5, 0.2, 0.3).normalized();
  Eigen::Matrix3d rot = q.toRotationMatrix();
  Eigen::Vector3d translation{0.6, 0.2, 0.43};
  raisim::Visuals* target = server_->addVisualSphere("target", 0.05);
  target->setPosition(translation);
  target->setColor(1.0, 0., 0., 0.92);
  pinocchio::SE3 desired_SE3(rot, translation);
  Eigen::VectorXd origin_q = pinocchio::neutral(controller_.getRobotModel());
  auto robotic_model = controller_.getRobotModel();
  auto robotic_data = controller_.getRobotData();
  int ee_frame_idx = controller_.getEndEffectorFrameId();
  ////////////////////////////////////////////////////////////////////////////////////////
  // Forward kinematics
  ////////////////////////////////////////////////////////////////////////////////////////
  pinocchio::framesForwardKinematics(robotic_model, robotic_data, origin_q);
  Eigen::Vector3d ee_pos_W = robotic_data.oMf[ee_frame_idx].translation();
  Eigen::Matrix3d ee_rot_W = robotic_data.oMf[ee_frame_idx].rotation();
  std::cout << "The EE position in world frame is: " << ee_pos_W.transpose()
            << std::endl;
  std::cout << "The EE rotation in world frame is: \n" << ee_rot_W << std::endl;
  auto testing_ee_pose = kinematic_chain_->forwardKinematic(origin_q);
  std::cout << "The testing EE position in world frame is : "
            << testing_ee_pose.block<3, 1>(0, 3).transpose() << std::endl;
  std::cout << "The testing EE rotation in world frame is: \n"
            << testing_ee_pose.block<3, 3>(0, 0) << std::endl;
  /////////////////////////////////////////////////////////////////////
  // Calculate jacobian
  /////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd jacobian_mat_B(6, robotic_model.nv),
      testing_jacobian_B(6, robotic_model.nv);
  auto start_time = std::chrono::system_clock::now();
  for (int i = 0; i < 10000; ++i)
    pinocchio::computeFrameJacobian(robotic_model, robotic_data, origin_q,
                                    ee_frame_idx, jacobian_mat_B);
  auto end_time = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);
  std::cout << "["
            << double(duration.count()) *
                   std::chrono::microseconds::period::num /
                   std::chrono::microseconds::period::den
            << "]" << "The Jacobian matrix is: \n"
            << jacobian_mat_B << std::endl;
  start_time = std::chrono::system_clock::now();
  for (int i = 0; i < 10000; ++i)
    kinematic_chain_->calculateBodyJacobian(origin_q, testing_jacobian_B);
  end_time = std::chrono::system_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                   start_time);
  std::cout << "["
            << double(duration.count()) *
                   std::chrono::microseconds::period::num /
                   std::chrono::microseconds::period::den
            << "]" << "The testing Jacobian matrix is: \n"
            << testing_jacobian_B << std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Inverse kinematic
  ///////////////////////////////////////////////////////////////////////////////////////////
  double error{0.0}, testing_error{0.0};
  Eigen::VectorXd testing_desired_theta(6);
  start_time = std::chrono::system_clock::now();
  Eigen::VectorXd desired_theta =
      controller_.solveInvKinematic(desired_SE3, origin_q, error);
  end_time = std::chrono::system_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                   start_time);
  std::cout
      << "["
      << double(duration.count()) * std::chrono::microseconds::period::num /
             std::chrono::microseconds::period::den
      << "]" << "(error: " << error << ")"
      << "Calculating inverse kinematics, the desired joint positions are: "
      << desired_theta.transpose() << std::endl;
  start_time = std::chrono::system_clock::now();
  kinematic_chain_->inverseKinematic(rot, translation, origin_q,
                                     testing_desired_theta, testing_error, 20);
  end_time = std::chrono::system_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                   start_time);
  std::cout << "["
            << double(duration.count()) *
                   std::chrono::microseconds::period::num /
                   std::chrono::microseconds::period::den
            << "]" << "(error: " << testing_error << ")"
            << "Testing calculating inverse kinematics, the testing desired "
               "joint positions are: "
            << testing_desired_theta.transpose() << std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////
  server_->launchServer();
  raisim::MSLEEP(1000);
  while (1) {
    robot_->setGeneralizedCoordinate(testing_desired_theta);
    // raisim::MSLEEP(4);
    // gf_.tail(num_joints_) =
    //     controller_.calculateTorque(desired_theta, gc_, gv_);
    // robot_->setGeneralizedForce(gf_);
    // if (server_) server_->lockVisualizationServerMutex();
    // world_->integrate();
    // if (server_) server_->unlockVisualizationServerMutex();
    // robot_->getState(gc_, gv_);
  }
}

}  // namespace demo
