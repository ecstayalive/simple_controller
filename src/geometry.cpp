#include "geometry.h"

#include <iostream>

void vecToSkewMat(const Eigen::Ref<const Eigen::Vector3d> &vec,
                  Eigen::Ref<Eigen::Matrix3d> skew_mat) {
  skew_mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
}
Eigen::Matrix3d vecToSkewMat(const Eigen::Ref<const Eigen::Vector3d> &vec) {
  Eigen::Matrix3d skew_mat;
  vecToSkewMat(vec, skew_mat);
  return skew_mat;
}

void skewMatToVec(const Eigen::Ref<const Eigen::Matrix3d> skew_mat,
                  Eigen::Ref<Eigen::Vector3d> vec) {
  vec << skew_mat(2, 1), skew_mat(0, 2), skew_mat(1, 0);
}
Eigen::Vector3d skewMatToVec(const Eigen::Ref<const Eigen::Matrix3d> skew_mat) {
  Eigen::Vector3d vec;
  skewMatToVec(skew_mat, vec);
  return vec;
}

void so3ToRot(const Eigen::Ref<const Eigen::Matrix3d> &skew_mat,
              const double &theta, Eigen::Ref<Eigen::Matrix3d> rot) {
  rot = Eigen::Matrix3d::Identity() + std::sin(theta) * skew_mat +
        (1 - std::cos(theta)) * skew_mat * skew_mat;
}
Eigen::Matrix3d so3ToRot(const Eigen::Ref<const Eigen::Matrix3d> &skew_mat,
                         const double &theta) {
  Eigen::Matrix3d rot;
  so3ToRot(skew_mat, theta, rot);
  return rot;
}

void angleAxisToRot(const Eigen::Ref<const Eigen::Vector3d> &axis_vec,
                    const double &theta, Eigen::Ref<Eigen::Matrix3d> rot) {
  Eigen::Matrix3d skew_mat;
  vecToSkewMat(axis_vec, skew_mat);
  rot = Eigen::Matrix3d::Identity() + std::sin(theta) * skew_mat +
        (1 - std::cos(theta)) * skew_mat * skew_mat;
}
Eigen::Matrix3d angleAxisToRot(
    const Eigen::Ref<const Eigen::Vector3d> &axis_vec, const double &theta) {
  Eigen::Matrix3d rot;
  angleAxisToRot(axis_vec, theta, rot);
  return rot;
}

void rotToAngleAxis(const Eigen::Ref<const Eigen::Matrix3d> &rot,
                    Eigen::Ref<Eigen::Vector3d> axis_vec, double &theta) {
  // theta belong to [0, pi]
  if ((rot - Eigen::Matrix3d::Identity()).norm() <= 1.0e-6) {
    // if rotation matrix is I
    theta = 0.0;
    // axis_vec is actually undefined
    // but for the sake of formal unity, assign it a unit vector
    axis_vec << 0., 0., 1.;
  } else if (std::abs(rot.trace() + 1.0) < 1e-5) {
    // if tr rot = -1, theta = pi
    theta = M_PI;
    Eigen::Vector3d bias_vec;
    bias_vec << 0, 0, 1;
    axis_vec = (bias_vec + rot.col(2)) / std::sqrt(2 * (1 + rot(2, 2)));
    // #ifdef DEBUG
    if (axis_vec.hasNaN() || std::abs(axis_vec.norm() - 1) > 0.05) {
      std::cout
          << "\033[1;31m rotToAngleAxis function error happens in 2th running "
             "case, axis_vec:\n"
          << axis_vec << "\nthe axis_vec's norm is: " << axis_vec.norm()
          << "\nrotation matrix:\n"
          << rot << "\nDenominator: " << std::sqrt(2 * (1 + rot(2, 2))) << "\n";
      throw std::runtime_error("rotToR3 function runtime error!!!\033[0m");
    }
    // #endif
  } else {
    theta = std::acos(std::max(std::min((rot.trace() - 1) * 0.5, 1.0), -1.0));
    Eigen::Matrix3d skew_mat =
        0.5 * (rot - rot.transpose()) / (std::sin(theta) + 1.0e-8);
    skewMatToVec(skew_mat, axis_vec);
    // #ifdef DEBUG
    if ((axis_vec * theta).hasNaN()) {
      std::cout
          << "\033[1;31m rotToAngleAxis function error happens in 3th running "
             "case, theta: "
          << theta << "\naxis_vec:\n"
          << axis_vec << "\naxis_vec's norm is: " << axis_vec.norm()
          << "\nrotation matrix:\n"
          << rot << "\ntrace of rotation matrix: " << rot.trace() << "\n";
      throw std::runtime_error("rotToR3 function runtime error!!!\033[0m");
    }
    // #endif
  }
}

void vecToRot(const Eigen::Ref<const Eigen::Vector3d> &vec,
              Eigen::Ref<Eigen::Matrix3d> rot) {
  double theta(vec.norm());
  if (std::abs(theta) < 1e-4) {
    rot = Eigen::Matrix3d::Identity();
  } else {
    Eigen::Vector3d axis_vec = vec / vec.norm();
    angleAxisToRot(axis_vec, theta, rot);
  }
}

Eigen::Matrix3d vecToRot(const Eigen::Ref<const Eigen::Vector3d> &vec) {
  Eigen::Matrix3d rot;
  vecToRot(vec, rot);
  return rot;
}

void rotToVec(const Eigen::Ref<const Eigen::Matrix3d> &rot,
              Eigen::Ref<Eigen::Vector3d> vec) {
  double theta;
  Eigen::Vector3d axis_vec;
  rotToAngleAxis(rot, axis_vec, theta);
  vec = theta * axis_vec;
}
Eigen::Vector3d rotToVec(const Eigen::Ref<const Eigen::Matrix3d> &rot) {
  Eigen::Vector3d vec;
  rotToVec(rot, vec);
  return vec;
}

void rotAndPosToTransMat(const Eigen::Ref<const Eigen::Matrix3d> &rot,
                         const Eigen::Ref<const Eigen::Vector3d> &pos,
                         Eigen::Ref<Eigen::Matrix4d> trans_mat) {
  trans_mat.block<3, 3>(0, 0) = rot;
  trans_mat.block<3, 1>(0, 3) = pos;
  trans_mat.block<1, 3>(3, 0).setZero();
  trans_mat(3, 3) = 1.0;
}
Eigen::Matrix4d rotAndPosToTransMat(
    const Eigen::Ref<const Eigen::Matrix3d> &rot,
    const Eigen::Ref<const Eigen::Vector3d> &pos) {
  Eigen::Matrix4d trans_mat;
  rotAndPosToTransMat(rot, pos, trans_mat);
  return trans_mat;
}

void transMatToRotAndPos(const Eigen::Ref<const Eigen::Matrix4d> &homo,
                         Eigen::Ref<Eigen::Matrix3d> rot,
                         Eigen::Ref<Eigen::Vector3d> pos) {
  rot = homo.block<3, 3>(0, 0);
  pos = homo.block<3, 1>(0, 3);
}

void twistToTwistHomo(const Eigen::Ref<const Eigen::VectorXd> &twist,
                      Eigen::Ref<Eigen::Matrix4d> twist_homo) {
  Eigen::Vector3d angular_vel = twist.head(3), linear_vel = twist.tail(3);
  twist_homo.setZero();
  twist_homo.block<3, 3>(0, 0) = vecToSkewMat(angular_vel);
  twist_homo.block<3, 1>(0, 3) = linear_vel;
}
Eigen::Matrix4d twistToTwistHomo(
    const Eigen::Ref<const Eigen::VectorXd> &twist) {
  Eigen::Matrix4d twist_homo;
  twistToTwistHomo(twist, twist_homo);
  return twist_homo;
}

void twistHomoToTwist(const Eigen::Ref<const Eigen::Matrix4d> &twist_homo,
                      Eigen::Ref<Eigen::VectorXd> twist) {
  Eigen::Matrix3d angular_skew_mat = twist_homo.block<3, 3>(0, 0);
  Eigen::Vector3d linear_vel = twist_homo.block<3, 1>(0, 3);
  twist.head(3) = skewMatToVec(angular_skew_mat);
  twist.tail(3) = linear_vel;
}
Eigen::VectorXd twistHomoToTwist(
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo) {
  Eigen::VectorXd twist(6);
  twistHomoToTwist(twist_homo, twist);
  return twist;
}

void adjointMat(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                Eigen::Ref<Eigen::MatrixXd> adjoint_mat) {
  Eigen::Matrix3d rot;
  Eigen::Vector3d pos;
  transMatToRotAndPos(trans_mat, rot, pos);
  adjoint_mat.setZero();
  adjoint_mat.block<3, 3>(0, 0) = rot;
  adjoint_mat.block<3, 3>(3, 0) = vecToSkewMat(pos) * rot;
  adjoint_mat.block<3, 3>(3, 3) = rot;
}

Eigen::MatrixXd adjointMat(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat) {
  Eigen::MatrixXd adjoint_mat(6, 6);
  adjointMat(trans_mat, adjoint_mat);
  return adjoint_mat;
}

void adjointMap(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                const Eigen::Ref<const Eigen::VectorXd> &twist,
                Eigen::Ref<Eigen::VectorXd> result) {
  // result = adjointMat(trans_mat) * twist;
  Eigen::Matrix4d trans_mat_inv = Eigen::Matrix4d::Identity();
  trans_mat_inv.block<3, 3>(0, 0) = trans_mat.block<3, 3>(0, 0).transpose();
  trans_mat_inv.block<3, 1>(0, 3) =
      -trans_mat_inv.block<3, 3>(0, 0) * trans_mat.block<3, 1>(0, 3);
  Eigen::Matrix4d twist_homo;
  twistToTwistHomo(twist, twist_homo);
  Eigen::Matrix4d result_homo = trans_mat * twist_homo * trans_mat_inv;
  twistHomoToTwist(result_homo, result);
}
Eigen::VectorXd adjointMap(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                           const Eigen::Ref<const Eigen::VectorXd> &twist) {
  Eigen::VectorXd result(6);
  adjointMap(trans_mat, twist, result);
  return result;
}

void twistHomoAdjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo,
    Eigen::Ref<Eigen::Matrix4d> result) {
  // result = adjointMat(trans_mat) * twist;
  Eigen::Matrix4d trans_mat_inv = Eigen::Matrix4d::Identity();
  trans_mat_inv.block<3, 3>(0, 0) = trans_mat.block<3, 3>(0, 0).transpose();
  trans_mat_inv.block<3, 1>(0, 3) =
      -trans_mat.block<3, 3>(0, 0).transpose() * trans_mat.block<3, 1>(0, 3);
  result = trans_mat * twist_homo * trans_mat_inv;
}

Eigen::Matrix4d twistHomoAdjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo) {
  Eigen::Matrix4d result;
  twistHomoAdjointMapToTwistHomo(trans_mat, twist_homo, result);
  return result;
}

void adjointMapToTwistHomo(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                           const Eigen::Ref<const Eigen::VectorXd> &twist,
                           Eigen::Ref<Eigen::Matrix4d> result) {
  // result = adjointMat(trans_mat) * twist;
  Eigen::Matrix4d trans_mat_inv = Eigen::Matrix4d::Identity();
  trans_mat_inv.block<3, 3>(0, 0) = trans_mat.block<3, 3>(0, 0).transpose();
  trans_mat_inv.block<3, 1>(0, 3) =
      -trans_mat.block<3, 3>(0, 0).transpose() * trans_mat.block<3, 1>(0, 3);
  Eigen::Matrix4d twist_homo;
  twistToTwistHomo(twist, twist_homo);
  result = trans_mat * twist_homo * trans_mat_inv;
}

Eigen::Matrix4d adjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::VectorXd> &twist) {
  Eigen::Matrix4d result;
  adjointMapToTwistHomo(trans_mat, twist, result);
  return result;
}

void twistToScrew(const Eigen::Ref<const Eigen::VectorXd> &twist,
                  Eigen::Ref<Eigen::VectorXd> screw_axis, double &d_theta) {
  Eigen::Vector3d angular_vel = twist.head(3), linear_vel = twist.head(3);
  if (std::abs(angular_vel.norm()) < 1.0e-6) {
    d_theta = linear_vel.norm();
  } else {
    d_theta = angular_vel.norm();
  }
  screw_axis = twist / d_theta;
}

void se3ToTransMat(const Eigen::Ref<const Eigen::Matrix4d> &screw_skew_mat,
                   const double &theta, Eigen::Ref<Eigen::Matrix4d> trans_mat) {
  Eigen::Matrix3d angular_skew_mat = screw_skew_mat.block<3, 3>(0, 0);
  Eigen::Vector3d vel_axis = screw_skew_mat.block<3, 1>(0, 3);
  trans_mat.setIdentity();
  so3ToRot(angular_skew_mat, theta, trans_mat.block<3, 3>(0, 0));
  trans_mat.block<3, 1>(0, 3) =
      (theta * Eigen::Matrix3d::Identity() +
       (1 - std::cos(theta)) * angular_skew_mat +
       (theta - std::sin(theta)) * angular_skew_mat * angular_skew_mat) *
      vel_axis;
}

Eigen::Matrix4d se3ToTransMat(
    const Eigen::Ref<const Eigen::Matrix4d> &screw_skew_mat,
    const double &theta) {
  Eigen::Matrix4d trans_mat;
  se3ToTransMat(screw_skew_mat, theta, trans_mat);
  return trans_mat;
}

void screwToTransMat(const Eigen::Ref<const Eigen::VectorXd> &screw_axis,
                     const double &theta,
                     Eigen::Ref<Eigen::Matrix4d> trans_mat) {
  Eigen::Vector3d angular_axis = screw_axis.head(3),
                  vel_axis = screw_axis.tail(3);
  Eigen::Matrix3d angular_skew_mat = vecToSkewMat(angular_axis);
  trans_mat.setIdentity();
  trans_mat.block<3, 3>(0, 0) = angleAxisToRot(angular_axis, theta);
  trans_mat.block<3, 1>(0, 3) =
      (theta * Eigen::Matrix3d::Identity() +
       (1 - std::cos(theta)) * angular_skew_mat +
       (theta - std::sin(theta)) * angular_skew_mat * angular_skew_mat) *
      vel_axis;
}
Eigen::Matrix4d screwToTransMat(
    const Eigen::Ref<const Eigen::VectorXd> &screw_axis, const double &theta) {
  Eigen::Matrix4d homo_mat;
  screwToTransMat(screw_axis, theta, homo_mat);
  return homo_mat;
}

void twistToTransMat(const Eigen::Ref<const Eigen::VectorXd> &twist,
                     Eigen::Ref<Eigen::Matrix4d> trans_mat) {
  double theta;
  Eigen::VectorXd screw_axis(6);
  twistToScrew(twist, screw_axis, theta);
  screwToTransMat(screw_axis, theta, trans_mat);
}
Eigen::Matrix4d twistToTransMat(
    const Eigen::Ref<const Eigen::VectorXd> &twist) {
  Eigen::Matrix4d trans_mat;
  twistToTransMat(twist, trans_mat);
  return trans_mat;
}

void transMatToScrew(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                     Eigen::Ref<Eigen::VectorXd> screw_axis, double &theta) {
  Eigen::Matrix3d rot;
  Eigen::Vector3d pos;
  transMatToRotAndPos(trans_mat, rot, pos);
  if ((rot - Eigen::Matrix3d::Identity()).norm() < 1.0e-6) {
    // if rotation matrix is identity
    // omega = 0
    screw_axis.head(3).setZero();
    theta = pos.norm();
    screw_axis.tail(3) = pos / std::max(pos.norm(), 1.0e-6);
    if (screw_axis.hasNaN()) {
      std::cout
          << "\033[1;31m transMatToScrew function error happens in 1th running "
             "case, theta: "
          << theta << "\nscrew_axis_vec:\n"
          << screw_axis << "\nscrew_axis's norm is: " << screw_axis.norm()
          << "\nrotation matrix:\n"
          << rot << "\ntranslation matrix: " << pos.transpose() << "\n";
      throw std::runtime_error(
          "transMatToScrew function runtime error!!!\033[0m");
    }
  } else {
    Eigen::Vector3d angular_axis;
    rotToAngleAxis(rot, angular_axis, theta);
    screw_axis.head(3) = angular_axis;
    Eigen::Matrix3d G_inv, angular_skew_mat;
    vecToSkewMat(angular_axis, angular_skew_mat);
    G_inv = Eigen::Matrix3d::Identity() / theta - 0.5 * angular_skew_mat +
            (1 / theta - 1 / (std::tan(theta / 2) * 2)) * angular_skew_mat *
                angular_skew_mat;
    screw_axis.tail(3) = G_inv * pos;
    if (screw_axis.hasNaN()) {
      std::cout
          << "\033[1;31m transMatToScrew function error happens in 2th running "
             "case, theta: "
          << theta << "\nscrew_axis_vec:\n"
          << screw_axis << "\nscrew_axis's norm is: " << screw_axis.norm()
          << "\nrotation matrix:\n"
          << rot << "\ntranslation matrix: " << pos.transpose() << "\n";
      throw std::runtime_error(
          "transMatToScrew function runtime error!!!\033[0m");
    }
  }
}

void transMatToTwist(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                     Eigen::Ref<Eigen::VectorXd> twist) {
  Eigen::VectorXd screw_axis(6);
  double theta;
  transMatToScrew(trans_mat, screw_axis, theta);
  twist = screw_axis * theta;
}
Eigen::VectorXd transMatToTwist(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat) {
  Eigen::VectorXd twist(6);
  transMatToTwist(trans_mat, twist);
  return twist;
}
