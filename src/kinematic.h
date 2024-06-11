#pragma once
#include <Eigen/Core>

struct TransMat {
  TransMat() : rot(Eigen::Matrix3d::Identity()), pos(Eigen::Vector3d::Zero()){};
  explicit TransMat(const Eigen::Ref<const Eigen::Matrix3d>& r,
                    const Eigen::Ref<const Eigen::Vector3d>& p)
      : rot(r), pos(p){};
  explicit TransMat(const Eigen::Ref<const Eigen::Matrix4d>& T)
      : rot(T.block<3, 3>(0, 0)), pos(T.block<3, 1>(0, 3)){};
  Eigen::Matrix3d rot;
  Eigen::Vector3d pos;
};

class KinematicChain {
 public:
  KinematicChain(std::vector<Eigen::Vector3d> joint_axis,
                 std::vector<Eigen::Matrix3d> joint_init_rot,
                 std::vector<Eigen::Vector3d> joint_init_pos,
                 Eigen::Vector3d end_effector_pos_P,
                 std::vector<double> joint_init_angles,
                 std::vector<Eigen::Vector2d> joint_angle_limits);
  ~KinematicChain() = default;
  TransMat forwardKinematic(Eigen::VectorXd joint_angles);
  TransMat forwardKinematicWithSpatialJacobian(
      Eigen::VectorXd joint_angles, Eigen::Ref<Eigen::MatrixXd> jacobian_B);
  TransMat forwardKinematicWithBodyJacobian(
      Eigen::VectorXd joint_angles, Eigen::Ref<Eigen::MatrixXd> jacobian_W);
  Eigen::MatrixXd calculateSpatialJacobian(Eigen::VectorXd joint_angles);
  Eigen::MatrixXd calculateBodyJacobian(Eigen::VectorXd joint_angles);
  bool inverseKinematic(Eigen::Matrix3d target_rot, Eigen::Vector3d target_pos,
                        Eigen::VectorXd init_guess,
                        Eigen::Ref<Eigen::VectorXd> joint_pos, double& error,
                        const unsigned int max_iterations = 10,
                        const double eps = 1.0e-3);

 protected:
  Eigen::Vector3d end_effector_pos_P_;
  std::vector<double> joint_init_angle_;
  Eigen::MatrixXd joint_angle_limits_;
  unsigned int num_joints_;
  ////////////////////////////////////////////////////////////////
  std::vector<Eigen::Matrix4d> screw_axis_homo_W_;
  std::vector<Eigen::Matrix4d> screw_axis_homo_P_;
  Eigen::Matrix4d M_, M_inv_;
  ///////////////////////////////////////////////////////////////
  Eigen::MatrixXd W_e_;
  Eigen::MatrixXd joints_identity_mat_;
  double lamba_{0.1};
};
