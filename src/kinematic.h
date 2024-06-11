#pragma once
#include <Eigen/Core>

class KinematicChain {
 public:
  explicit KinematicChain(
      const std::vector<Eigen::Vector3d>& joint_axis,
      const std::vector<Eigen::Matrix3d>& joint_init_rot,
      const std::vector<Eigen::Vector3d>& joint_init_pos,
      const Eigen::Ref<const Eigen::Vector3d>& end_effector_pos_P,
      const Eigen::Ref<const Eigen::VectorXd>& joint_init_angles,
      const std::vector<Eigen::Vector2d>& joint_angle_limits);
  ~KinematicChain() = default;
  /**
   * @brief
   * @param joint_angles
   * @param target_pose
   */
  void forwardKinematic(const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
                        Eigen::Ref<Eigen::Matrix4d> target_pose);
  /**
   * @brief
   * @param joint_angles
   * @return Eigen::Matrix4d
   */
  Eigen::Matrix4d forwardKinematic(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles);
  /**
   * @brief
   * @param joint_angles
   * @param target_pose
   * @param spatial_jacobian
   */
  void forwardKinematicWithSpatialJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::Matrix4d> target_pose,
      Eigen::Ref<Eigen::MatrixXd> spatial_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @param spatial_jacobian
   * @return Eigen::Matrix4d Target pose
   */
  Eigen::Matrix4d forwardKinematicWithSpatialJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::MatrixXd> spatial_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @param body_jacobian
   * @return Eigen::Matrix4d
   */
  void forwardKinematicWithBodyJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::Matrix4d> target_pose,
      Eigen::Ref<Eigen::MatrixXd> body_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @param body_jacobian
   * @return Eigen::Matrix4d
   */
  Eigen::Matrix4d forwardKinematicWithBodyJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::MatrixXd> body_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @param spatial_jacobian
   */
  void calculateSpatialJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::MatrixXd> spatial_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd calculateSpatialJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles);
  /**
   * @brief
   * @param joint_angles
   * @param body_jacobian
   */
  void calculateBodyJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles,
      Eigen::Ref<Eigen::MatrixXd> body_jacobian);
  /**
   * @brief
   * @param joint_angles
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd calculateBodyJacobian(
      const Eigen::Ref<const Eigen::VectorXd>& joint_angles);
  /**
   * @brief
   * @param target_rot
   * @param target_pos
   * @param init_guess
   * @param joint_pos
   * @param error
   * @param max_iterations
   * @param eps
   * @return true
   * @return false
   */
  bool inverseKinematic(const Eigen::Ref<const Eigen::Matrix3d>& target_rot,
                        const Eigen::Ref<const Eigen::Vector3d>& target_pos,
                        const Eigen::Ref<const Eigen::VectorXd>& init_guess,
                        Eigen::Ref<Eigen::VectorXd> joint_pos, double& error,
                        const unsigned int max_iterations = 10,
                        const double eps = 1.0e-3);
  /**
   * @brief
   * @param target_pose
   * @param init_guess
   * @param joint_pos
   * @param error
   * @param max_iterations
   * @param eps
   * @return true
   * @return false
   */
  bool inverseKinematic(const Eigen::Ref<const Eigen::Matrix4d>& target_pose,
                        const Eigen::Ref<const Eigen::VectorXd>& init_guess,
                        Eigen::Ref<Eigen::VectorXd> joint_pos, double& error,
                        const unsigned int max_iterations = 10,
                        const double eps = 1.0e-3);
  /**
   * @brief Set the Inverse Kinematic Solver Params object
   * @param lamba
   * @param error_wight
   */
  void setInvKinematicSolverParams(
      const double& lamba,
      const Eigen::Ref<const Eigen::VectorXd>& error_wight) {
    lamba_ = lamba;
    error_wight_ = error_wight;
  }

 protected:
  // kinematic chain params
  unsigned int num_joints_;
  Eigen::VectorXd joint_init_angles_;
  Eigen::MatrixXd joint_angles_limits_;
  std::vector<Eigen::Matrix4d> screw_axis_homo_W_;
  std::vector<Eigen::Matrix4d> screw_axis_homo_B_;
  Eigen::Matrix4d init_trans_mat_, init_trans_mat_inv_;
  // inverse kinematic solver params
  Eigen::MatrixXd kJointIdentityMat_;
  Eigen::VectorXd error_wight_;  // TODO: Unused
  double lamba_{0.1};
};
