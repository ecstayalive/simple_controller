#include "controller.h"

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/urdf.hpp"

namespace demo {

void Controller::init(const std::string& urdf_path) {
  pinocchio::urdf::buildModel(urdf_path, model_);
  pinocchio::Data data(model_);
  end_effector_frame_id_ =
      model_.getFrameId("endEffector", pinocchio::FIXED_JOINT);
  data_ = data;
  identity_mat_.setIdentity();
  kp_.setIdentity(model_.nv, model_.nv);
  kd_.setIdentity(model_.nv, model_.nv);
  kp_ = kp_ * 30;
  kd_ = kd_ * 1.0;
}

Eigen::VectorXd Controller::solveInvKinematic(pinocchio::SE3 desired_pose,
                                              Eigen::VectorXd initial_guess,
                                              double& error) {
  Eigen::MatrixXd jacobian_mat_B(
      6, model_.nv);  ///< jacobian matrix under body frame
  Eigen::MatrixXd jacobian_mat_E(
      6, model_.nv);  ///< jacobian matrix under error task
  Eigen::MatrixXd jaco_log6(6, model_.nv);
  Eigen::VectorXd v(model_.nv);
  Eigen::VectorXd error_vec_B(6);
  bool success{false};
  Eigen::VectorXd theta = initial_guess;
  int iteration = 0;
  while (iteration < max_iterations_) {
    // pinocchio::forwardKinematics(model_, data_, theta);
    // // Get the end-effector pose under body frame
    // pinocchio::SE3 pose_B = data_.oMi[end_joint_id_].actInv(desired_pose);
    pinocchio::framesForwardKinematics(model_, data_, theta);
    pinocchio::SE3 pose_B =
        data_.oMf[end_effector_frame_id_].actInv(desired_pose);
    error_vec_B = pinocchio::log6(pose_B).toVector();
    if (error_vec_B.norm() < eps_) {
      success = true;
      break;
    }
    // Calculate body jacobian
    // pinocchio::computeJointJacobian(model_, data_, theta, end_joint_id_,
    //                                 jacobian_mat_B);
    pinocchio::computeFrameJacobian(model_, data_, theta,
                                    end_effector_frame_id_, jacobian_mat_B);
    // pinocchio::Jlog6(pose_B.inverse(), jaco_log6);
    // jacobian_mat_E = -jaco_log6 * jacobian_mat_B;
    // pinocchio::Data::Matrix6 JJt;
    // JJt.noalias() = jacobian_mat_E * jacobian_mat_E.transpose();
    // v.noalias() = -jacobian_mat_E.transpose() *
    // JJt.ldlt().solve(error_vec_B);
    pinocchio::Jlog6(pose_B, jaco_log6);
    jacobian_mat_E = -jaco_log6 * jacobian_mat_B;
    // use pseudo inverse matrix
    v = -jacobian_mat_E.transpose() *
        (jacobian_mat_E * jacobian_mat_E.transpose() + damp_ * identity_mat_)
            .inverse() *
        error_vec_B;
    theta = pinocchio::integrate(model_, theta, v * dt_);
    ++iteration;
  }
  error = error_vec_B.norm();
  return theta;
}

Eigen::VectorXd Controller::calculateTorque(Eigen::VectorXd desired_theta,
                                            Eigen::VectorXd q,
                                            Eigen::VectorXd dq) {
  Eigen::VectorXd theta_error = desired_theta - q;
  Eigen::VectorXd d_theta_error = -dq;
  // feedback control
  Eigen::VectorXd desired_theta_acc = kp_ * theta_error + kd_ * d_theta_error;
  Eigen::VectorXd joint_torque;
  // feedforward control
  joint_torque = pinocchio::rnea(model_, data_, q, dq, desired_theta_acc);
  return joint_torque;
}

}  // namespace demo
