#pragma once

#include <Eigen/Core>

#include "pinocchio/multibody/data.hpp"
#include "pinocchio/multibody/fwd.hpp"

namespace demo {

class Controller {
 public:
  void init(const std::string &urdf_path);

  /**
   * @brief Use task error Jacobian matrix to solve the inverse kinematic
   * problem
   * @version 2.0
   * @param desired_pose
   * @param initial_guess
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd solveInvKinematic(Eigen::Matrix4d desired_pose_mat,
                                    Eigen::VectorXd initial_guess) {
    pinocchio::SE3 desired_pose(desired_pose_mat.block<3, 3>(0, 0),
                                desired_pose_mat.block<3, 1>(0, 3));
    return solveInvKinematic(desired_pose, initial_guess);
  }

  /**
   * @brief Use task error Jacobian matrix to solve the inverse kinematic
   * problem
   * @version 2.0
   * @param desired_pose
   * @param initial_guess
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd solveInvKinematic(pinocchio::SE3 desired_pose,
                                    Eigen::VectorXd initial_guess);
  /**
   * @brief Use body pose Jacobian matrix to solve the inverse kinematic problem
   * @details Compared with using the error task Jacobian matrix, this method
   * has lower applicability and is accompanied by instability in the
   * intermediate process.
   * @version 1.0
   * @param desired_pose
   * @param initial_guess
   * @return Eigen::VectorXd The result
   */
  Eigen::VectorXd solveInvKinematic1(pinocchio::SE3 desired_pose,
                                     Eigen::VectorXd initial_guess);

  /**
   * @brief
   *
   * @param desired_theta
   * @param q
   * @param dq
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd calculateTorque(Eigen::VectorXd desired_theta,
                                  Eigen::VectorXd q, Eigen::VectorXd dq);

  pinocchio::Model getRobotModel() const { return model_; }

  void setMaxIterations(int max_iterations) {
    max_iterations_ = max_iterations;
  }
  void setEndJointId(int joint_id) { end_joint_id_ = joint_id; }

 protected:
  pinocchio::Model model_;
  pinocchio::Data data_;
  int end_joint_id_ = 6;
  int end_effector_frame_id_;
  ////////////////////////////////
  // inverse kinematics
  ////////////////////////////////
  double eps_ = 1e-4;  ///< the twist error
  int max_iterations_ = 1000;
  double dt_ = 0.1;
  double damp_ = 1e-3;
  Eigen::Matrix<double, 6, 6> identity_mat_;
  ////////////////////////////////
  // controller
  ////////////////////////////////
  Eigen::MatrixXd kp_;
  Eigen::MatrixXd kd_;
};

}  // namespace demo
