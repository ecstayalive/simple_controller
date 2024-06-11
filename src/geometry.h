#pragma once
#include <Eigen/Core>

/**
 * @brief Transform the vector to its skew symmetric matrix
 *
 * @param vec
 * @param skew_mat
 */
void vecToSkewMat(const Eigen::Ref<const Eigen::Vector3d> &vec,
                  Eigen::Ref<Eigen::Matrix3d> skew_mat);
[[nodiscard]] Eigen::Matrix3d vecToSkewMat(
    const Eigen::Ref<const Eigen::Vector3d> &vec);

void skewMatToVec(const Eigen::Ref<const Eigen::Matrix3d> skew_mat,
                  Eigen::Ref<Eigen::Vector3d> vec);
[[nodiscard]] Eigen::Vector3d skewMatToVec(
    const Eigen::Ref<const Eigen::Matrix3d> skew_mat);

/**
 * @brief Mapping of so(3) space to rotation matrix SO(3)
 *
 * @param skew_mat the skew symmetric matrix which is converted from
              the unit vector
 * @param theta the rotation angle which belongs to (-pi, pi]
 * @param rot rotation matrix
 */
void so3ToRot(const Eigen::Ref<const Eigen::Matrix3d> &skew_mat,
              const double &theta, Eigen::Ref<Eigen::Matrix3d> rot);
[[nodiscard]] Eigen::Matrix3d so3ToRot(
    const Eigen::Ref<const Eigen::Matrix3d> &skew_mat, const double &theta);

/**
 * @brief Mapping of r(3) space to rotation matrix SO(3)
 *
 * @param axis_vec the rotation axis, is presented as a unit vector
 * @param theta the rotation angle which belongs to (-pi, pi]
 * @param rot rotation matrix
 */
void angleAxisToRot(const Eigen::Ref<const Eigen::Vector3d> &axis_vec,
                    const double &theta, Eigen::Ref<Eigen::Matrix3d> rot);
[[nodiscard]] Eigen::Matrix3d angleAxisToRot(
    const Eigen::Ref<const Eigen::Vector3d> &axis_vec, const double &theta);

/**
 * @brief Mapping of rotation matrix SO(3) to r(3) space
 *
 * @param rot rotation matrix
 * @param axis_vec the rotation axis, is presented as a unit vector
 * @param theta the rotation angle which belongs to [0, pi]
 */
void rotToAngleAxis(const Eigen::Ref<const Eigen::Matrix3d> &rot,
                    Eigen::Ref<Eigen::Vector3d> axis_vec, double &theta);

/**
 * @brief Mapping of R(3) space to rotation matrix SO(3)
 *
 * @param vec the L2 norm is the rotation angle and
              the direction is the rotation axis
 * @param rot rotation matrix
 */
void vecToRot(const Eigen::Ref<const Eigen::Vector3d> &vec,
              Eigen::Ref<Eigen::Matrix3d> rot);
[[nodiscard]] Eigen::Matrix3d vecToRot(
    const Eigen::Ref<const Eigen::Vector3d> &vec);

/**
 * @brief Mapping of SO(3) space to R(3)
 *
 * @param rot rotation matrix
 * @param vec the L2 norm is the rotation angle and
              the direction is the rotation axis
 */
void rotToVec(const Eigen::Ref<const Eigen::Matrix3d> &rot,
              Eigen::Ref<Eigen::Vector3d> vec);
[[nodiscard]] Eigen::Vector3d rotToVec(
    const Eigen::Ref<const Eigen::Matrix3d> &rot);

/**
 * @brief Transform rotation matrix and position vector
          to homogeneous transformation matrix.
 * @param rot rotation matrix
 * @param pos position vector
 * @param trans_mat homogeneous transformation matrix
 */
void rotAndPosToTransMat(const Eigen::Ref<const Eigen::Matrix3d> &rot,
                         const Eigen::Ref<const Eigen::Vector3d> &pos,
                         Eigen::Ref<Eigen::Matrix4d> trans_mat);
[[nodiscard]] Eigen::Matrix4d rotAndPosToTransMat(
    const Eigen::Ref<const Eigen::Matrix3d> &rot,
    const Eigen::Ref<const Eigen::Vector3d> &pos);
/**
 * @brief Transform homogeneous transformation matrix to
          rotation matrix and position vector.
 * @param trans_mat
 * @param rot
 * @param pos
 */
void transMatToRotAndPos(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                         Eigen::Ref<Eigen::Matrix3d> rot,
                         Eigen::Ref<Eigen::Vector3d> pos);

/**
 * @brief Get 4x4 matrix representation of twist
 *
 * @param twist 6x1 vector
 * @param twist_homo_mat 4x4 matrix
 */
void twistToTwistHomo(const Eigen::Ref<const Eigen::VectorXd> &twist,
                      Eigen::Ref<Eigen::Matrix4d> twist_homo);
[[nodiscard]] Eigen::Matrix4d twistToTwistHomo(
    const Eigen::Ref<const Eigen::VectorXd> &twist);

/**
 * @brief Get twist from its homogeneous matrix representation
 *
 * @param twist_homo 4x4 matrix
 * @param twist 6x1 vector
 */
void twistHomoToTwist(const Eigen::Ref<const Eigen::Matrix4d> &twist_homo,
                      Eigen::Ref<Eigen::VectorXd> twist);
[[nodiscard]] Eigen::VectorXd twistHomoToTwist(
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo);

/**
 * @brief Get homogeneous transformation's adjoint presentation
 *
 * @param trans_mat
 * @param adjoint_mat
 */
void adjointMat(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                Eigen::Ref<Eigen::MatrixXd> adjoint_mat);
[[nodiscard]] Eigen::MatrixXd adjointMat(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat);

/**
 * @brief Change representation coordinates of one twist
 *
 * @param trans_mat homogeneous transformation matrix
 * @param twist
 * @param result
 */
void adjointMap(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                const Eigen::Ref<const Eigen::VectorXd> &twist1,
                Eigen::Ref<Eigen::VectorXd> result);
[[nodiscard]] Eigen::VectorXd adjointMap(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::VectorXd> &twist);

/**
 * @brief Change representation coordinates of one twist
 *
 * @param trans_mat homogeneous transformation matrix
 * @param twist_homo
 * @param result
 */
void twistHomoAdjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo,
    Eigen::Ref<Eigen::Matrix4d> result);
[[nodiscard]] Eigen::Matrix4d twistHomoAdjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::Matrix4d> &twist_homo);

/**
 * @brief Change representation coordinates of one twist
 *
 * @param trans_mat homogeneous transformation matrix
 * @param twist
 * @param result
 */
void adjointMapToTwistHomo(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                           const Eigen::Ref<const Eigen::VectorXd> &twist,
                           Eigen::Ref<Eigen::Matrix4d> result);
[[nodiscard]] Eigen::Matrix4d adjointMapToTwistHomo(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
    const Eigen::Ref<const Eigen::VectorXd> &twist);

/**
 * @brief Use screw to represent twist
 * @details This also can be used to convert one screw vector to its axis and
 * norm
 * @param twist
 * @param screw_axis 6x1 vector
 * @param d_theta derivation of theta
 */
void twistToScrew(const Eigen::Ref<const Eigen::VectorXd> &twist,
                  Eigen::Ref<Eigen::VectorXd> screw_axis, double &d_theta);

/**
 * @brief Convert screw to homogeneous transformation(R6 -> se3 -> SE3)
 *
 * @param screw_skew_mat screw axis (6x1 vector)
 * @param theta
 * @param trans_mat homogeneous transformation
 */
void se3ToTransMat(const Eigen::Ref<const Eigen::Matrix4d> &screw_skew_mat,
                   const double &theta, Eigen::Ref<Eigen::Matrix4d> trans_mat);
[[nodiscard]] Eigen::Matrix4d se3ToTransMat(
    const Eigen::Ref<const Eigen::Matrix4d> &screw_skew_mat,
    const double &theta);

/**
 * @brief Convert screw to homogeneous transformation(R6 -> se3 -> SE3)
 *
 * @param screw_axis screw axis (6x1 vector)
 * @param theta
 * @param trans_mat homogeneous transformation
 */
void screwToTransMat(const Eigen::Ref<const Eigen::VectorXd> &screw_axis,
                     const double &theta,
                     Eigen::Ref<Eigen::Matrix4d> trans_mat);
[[nodiscard]] Eigen::Matrix4d screwToTransMat(
    const Eigen::Ref<const Eigen::VectorXd> &screw_axis, const double &theta);

/**
 * @brief  Convert screw vector to homogeneous transformation(R6 -> se3 -> SE3)
 *
 * @param screw_vec screw axis * theta
 * @param trans_mat homogeneous transformation
 */
void twistToTransMat(const Eigen::Ref<const Eigen::VectorXd> &twist,
                     Eigen::Ref<Eigen::Matrix4d> trans_mat);
[[nodiscard]] Eigen::Matrix4d twistToTransMat(
    const Eigen::Ref<const Eigen::VectorXd> &twist);

/**
 * @brief Convert homogeneous transformation to screw (SE3 -> se3 -> R6)
 *
 * @param trans_mat
 * @param screw_axis
 * @param theta
 */
void transMatToScrew(const Eigen::Ref<const Eigen::Matrix4d> &homo_mat,
                     Eigen::Ref<Eigen::VectorXd> screw_axis, double &theta);

/**
 * @brief Convert homogeneous transformation to screw vector (SE3 -> se3 -> R6)
 *
 * @param trans_mat
 * @param screw_vec
 */
void transMatToTwist(const Eigen::Ref<const Eigen::Matrix4d> &trans_mat,
                     Eigen::Ref<Eigen::VectorXd> twist);
[[nodiscard]] Eigen::VectorXd transMatToTwist(
    const Eigen::Ref<const Eigen::Matrix4d> &trans_mat);
