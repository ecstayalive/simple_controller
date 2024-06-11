#include "kinematic.h"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>

#include "geometry.h"

KinematicChain::KinematicChain(std::vector<Eigen::Vector3d> joint_axis,
                               std::vector<Eigen::Matrix3d> joint_init_rot,
                               std::vector<Eigen::Vector3d> joint_init_pos,
                               Eigen::Vector3d end_effector_pos_P,
                               std::vector<double> joint_init_angle,
                               std::vector<Eigen::Vector2d> joint_angle_limits)
    : end_effector_pos_P_(end_effector_pos_P),
      joint_init_angle_(joint_init_angle) {
  num_joints_ = joint_axis.size();
  assert(num_joints_ == joint_init_pos.size());
  assert(num_joints_ == joint_init_rot.size());
  assert(num_joints_ == joint_init_angle.size());
  assert(num_joints_ == joint_angle_limits.size());
  std::cout << "The number of joints is: " << num_joints_ << std::endl;
  joint_angle_limits_.setZero(num_joints_, 2);
  M_.setIdentity();
  for (int i = 0; i < num_joints_; ++i) {
    joint_angle_limits_.row(i) = joint_angle_limits[i];
    Eigen::Matrix4d trans_mat = Eigen::Matrix4d::Identity();
    Eigen::VectorXd axis(6);
    Eigen::Matrix4d screw_axis_homo_W;
    axis.setZero();
    axis.head(3) = joint_axis[i];
    trans_mat.block<3, 3>(0, 0) = joint_init_rot[i];
    trans_mat.block<3, 1>(0, 3) = joint_init_pos[i];
    M_ = M_ * trans_mat;
    adjointMapToTwistHomo(M_, axis, screw_axis_homo_W);
    screw_axis_homo_W_.emplace_back(screw_axis_homo_W);
  }
  M_.block<3, 1>(0, 3) =
      M_.block<3, 1>(0, 3) + M_.block<3, 3>(0, 0) * end_effector_pos_P;
  M_inv_.setIdentity();
  M_inv_.block<3, 3>(0, 0) = M_.block<3, 3>(0, 0).transpose();
  M_inv_.block<3, 1>(0, 3) = -M_inv_.block<3, 3>(0, 0) * M_.block<3, 1>(0, 3);
  for (int i = 0; i < num_joints_; ++i) {
    screw_axis_homo_P_.emplace_back(
        twistHomoAdjointMapToTwistHomo(M_inv_, screw_axis_homo_W_[i]));
  }

  W_e_.setIdentity(num_joints_, num_joints_);
  joints_identity_mat_.setIdentity(num_joints_, num_joints_);
}

TransMat KinematicChain::forwardKinematic(Eigen::VectorXd joint_angles) {
  assert(joint_angles.size() == num_joints_);
  Eigen::Matrix4d trans_mat = Eigen::Matrix4d::Identity();
  for (int i = 0; i < num_joints_; ++i) {
    trans_mat =
        trans_mat * se3ToTransMat(screw_axis_homo_W_[i], joint_angles[i]);
  }
  TransMat pose(trans_mat * M_);
  return pose;
}

TransMat KinematicChain::forwardKinematicWithSpatialJacobian(
    Eigen::VectorXd joint_angles,
    Eigen::Ref<Eigen::MatrixXd> spatial_jacobian) {
  assert(joint_angles.size() == num_joints_);
  Eigen::Matrix4d trans_mat = Eigen::Matrix4d::Identity();
  for (int i = 0; i < num_joints_; ++i) {
    spatial_jacobian.col(i) = twistHomoToTwist(
        twistHomoAdjointMapToTwistHomo(trans_mat, screw_axis_homo_W_[i]));
    trans_mat =
        trans_mat * se3ToTransMat(screw_axis_homo_W_[i], joint_angles[i]);
  }
  TransMat pose(trans_mat * M_);
  return pose;
}

TransMat KinematicChain::forwardKinematicWithBodyJacobian(
    Eigen::VectorXd joint_angles, Eigen::Ref<Eigen::MatrixXd> body_jacobian) {
  assert(joint_angles.size() == num_joints_);
  Eigen::Matrix4d trans_mat = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d trans_mat_inv = Eigen::Matrix4d::Identity();
  for (int i = num_joints_ - 1; i >= 0; --i) {
    trans_mat_inv.block<3, 3>(0, 0) = trans_mat.block<3, 3>(0, 0).transpose();
    trans_mat_inv.block<3, 1>(0, 3) =
        -trans_mat_inv.block<3, 3>(0, 0) * trans_mat.block<3, 1>(0, 3);
    body_jacobian.col(i) = twistHomoToTwist(
        twistHomoAdjointMapToTwistHomo(trans_mat_inv, screw_axis_homo_P_[i]));
    trans_mat =
        se3ToTransMat(screw_axis_homo_P_[i], joint_angles[i]) * trans_mat;
  }
  TransMat pose(M_ * trans_mat);
  return pose;
}

Eigen::MatrixXd KinematicChain::calculateSpatialJacobian(
    Eigen::VectorXd joint_angles) {
  assert(joint_angles.size() == num_joints_);
  Eigen::MatrixXd jacobian_W(6, num_joints_);
  Eigen::Matrix4d jacobian_trans_mat_W = Eigen::Matrix4d::Identity();
  for (int i = 0; i < num_joints_; ++i) {
    jacobian_W.col(i) = twistHomoToTwist(twistHomoAdjointMapToTwistHomo(
        jacobian_trans_mat_W, screw_axis_homo_W_[i]));
    jacobian_trans_mat_W =
        jacobian_trans_mat_W *
        se3ToTransMat(screw_axis_homo_W_[i], joint_angles[i]);
  }
  return jacobian_W;
}
Eigen::MatrixXd KinematicChain::calculateBodyJacobian(
    Eigen::VectorXd joint_angles) {
  assert(joint_angles.size() == num_joints_);
  Eigen::MatrixXd jacobian_P(6, num_joints_);
  Eigen::Matrix4d jacobian_trans_mat_P = Eigen::Matrix4d::Identity();
  for (int i = num_joints_ - 1; i >= 0; --i) {
    jacobian_P.col(i) = twistHomoToTwist(twistHomoAdjointMapToTwistHomo(
        jacobian_trans_mat_P, screw_axis_homo_P_[i]));
    jacobian_trans_mat_P =
        jacobian_trans_mat_P *
        se3ToTransMat(-screw_axis_homo_P_[i], joint_angles[i]);
  }
  return jacobian_P;
}

bool KinematicChain::inverseKinematic(
    Eigen::Matrix3d target_rot, Eigen::Vector3d target_pos,
    Eigen::VectorXd init_guess, Eigen::Ref<Eigen::VectorXd> joint_pos,
    double& error, const unsigned int max_iterations, const double eps) {
  assert(init_guess.size() == num_joints_);
  bool success = false;
  Eigen::MatrixXd jacobian_P(6, num_joints_);
  Eigen::Matrix4d target_pose = Eigen::Matrix4d::Identity();
  target_pose.block<3, 3>(0, 0) = target_rot;
  target_pose.block<3, 1>(0, 3) = target_pos;
  TransMat pose;
  Eigen::Matrix4d pose_inv = Eigen::Matrix4d::Identity();
  Eigen::VectorXd error_twist(6);
  joint_pos = init_guess;
  unsigned int iter = 0;
  while (iter < max_iterations) {
    pose = forwardKinematicWithBodyJacobian(joint_pos, jacobian_P);
    pose_inv.block<3, 3>(0, 0) = pose.rot.transpose();
    pose_inv.block<3, 1>(0, 3) = -pose_inv.block<3, 3>(0, 0) * pose.pos;
    transMatToTwist(pose_inv * target_pose, error_twist);
    if (error_twist.norm() < eps) {
      success = true;
      break;
    }
    Eigen::MatrixXd G_x = jacobian_P.transpose() * error_twist;
    Eigen::MatrixXd A_x =
        jacobian_P.transpose() * jacobian_P +
        lamba_ * 0.5 * error_twist.squaredNorm() * joints_identity_mat_;
    joint_pos = (joint_pos + A_x.llt().solve(G_x))
                    .cwiseMax(joint_angle_limits_.col(0))
                    .cwiseMin(joint_angle_limits_.col(1));
    ++iter;
  }
  error = error_twist.norm();
  return success;
}
