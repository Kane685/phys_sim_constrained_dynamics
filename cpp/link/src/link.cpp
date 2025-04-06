#include "link/include/link.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

Link::Link(const LinkType type)
    : type_(type), index_(0), mass_(0), inv_mass_(0),
    inertia_(Vector3r::Zero()), inv_inertia_(Vector3r::Zero()),
    vertex_num_(0), vertices_(3, 0), element_num_(0), elements_(3, 0),
    translation_(Vector3r::Zero()),
    rotation_(Matrix3r::Zero()),
    linear_velocity_(Vector3r::Zero()),
    angular_velocity_(Vector3r::Zero()),
    force_(Vector3r::Zero()),
    torque_(Vector3r::Zero()) {}

void Link::ApplyForce(const Vector3r& attach_point, const Vector3r& force) {
    force_ += force;
    torque_ += (rotation_ * attach_point).cross(force);
}

void Link::Step(const real time_step, const Options& opt) {
    const std::string error_location = "link::Link::Step";
    CheckCondition(time_step > 0, error_location,
        "Time step must be positive.");

    const real h = time_step;

    // Update velocity.
    linear_velocity_ += inv_mass_ * force_ * h;
    angular_velocity_ += ToInvWorldInertia() * (torque_ -
        angular_velocity_.cross(ToWorldInertia() * angular_velocity_)) * h;

    // Update position.
    translation_ += linear_velocity_ * h;
    rotation_ = ToRotationMatrix(angular_velocity_ * h) * rotation_;
}

const Matrix3r Link::ToWorldInertia() const {
    ////////////////////////////////////////////////////////////////////////////
    // Task 2.1 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Implement the inertia tensor in the world frame with the rotation matrix
    // rotation_.
    // We have ensured the link's inertia tensor in the body frame is a diagonal
    // matrix whose diagonal entries are inertia_, a 3D vector.
    // Refer to the lecture slides or the supplementary reading materials to
    // figure out how to convert the inertia tensor to the world frame.
    //
    // TODO.

    Matrix3r I_b = inertia_.asDiagonal();
    Matrix3r I_w = rotation_ * I_b * rotation_.transpose();

    return I_w;
}

const Matrix3r Link::ToInvWorldInertia() const {
    return rotation_ * inv_inertia_.asDiagonal() * rotation_.transpose();
}

void Link::set_R(const Matrix3r& R) {
    const Matrix3r I = Matrix3r::Identity();
    CheckCondition(IsRotationMatrix(R), "link::Link::set_R",
        "The rotation matrix is not orthogonal.");
    rotation_ = R;
}

void Link::set_index(const integer index) {
    CheckCondition(index >= 0, "link::Link::set_index", "Invalid index.");
    index_ = index;
}

void Link::set_mass_and_inv_mass(const real mass) {
    CheckCondition(mass > 0, "link::Link::set_mass_and_inv_mass",
        "Mass must be strictly positive.");
    mass_ = mass;
    inv_mass_ = ToReal(1) / mass;
}

void Link::set_inertia_and_inv_inertia(const Vector3r& inertia) {
    CheckCondition(inertia.minCoeff() > 0,
        "link::Link::set_inertia_and_inv_inertia",
        "Inertia must be strictly positive.");
    inertia_ = inertia;
    inv_inertia_ = inertia.cwiseInverse();
}

void Link::set_vertices(const Matrix3Xr& vertices) {
    vertex_num_ = static_cast<integer>(vertices.cols());
    CheckCondition(vertex_num_ >= 3, "link::Link::set_vertices",
        "Need at least 3 vertices.");
    vertices_ = vertices;
}

void Link::set_elements(const Matrix3Xi& elements) {
    element_num_ = static_cast<integer>(elements.cols());
    CheckCondition(element_num_ >= 1, "link::Link::set_elements",
        "Need at least 1 triangle.");
    elements_ = elements;
}

const Eigen::Matrix<real, 3, 6> Link::ComputePointJacobian(
    const Vector3r& point) const {
    ////////////////////////////////////////////////////////////////////////////
    // Task 1.1 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // The input point is the 3D coordinates of a point in the body frame of the
    // rigid link. This function returns a 3 x 6 matrix J so that
    //
    // J (v)
    //   (w)
    //
    // computes the velocity of this point in the world frame. Here, v and w
    // are the linear and angular velocities of the link.
    //
    // TODO.

    // define J
    Eigen::Matrix<real, 3, 6> J;
    // J_v = I
    J.block<3,3>(0,0) = Matrix3r::Identity();
    // J_w = -[Rb]
    Vector3r Rb = rotation_ * point;
    // combine J_v and J_w
    J(0,3) =     0;    J(0,4) = -Rb(2); J(0,5) =  Rb(1);
    J(1,3) =  Rb(2);   J(1,4) =     0;  J(1,5) = -Rb(0);
    J(2,3) = -Rb(1);   J(2,4) =  Rb(0); J(2,5) =     0;

    return J;
}

const Eigen::Matrix<real, 3, 6> Link::ComputeVectorJacobian(
    const Vector3r& vector) const {

    // Consider phi = R p.
    // phi_dot = J * (v, w) = w x (Rp) = -[Rp] w.
    // Therefore,
    // J = (0,  [-Rp]).

    Eigen::Matrix<real, 3, 6> J; J.setZero();
    J.rightCols(3) = -ToCrossProductMatrix(ToWorldVector(vector));
    return J;
}

const std::string ToString(const LinkType type) {
    switch (type) {
        case LinkType::kBox:
            return "box";
        case LinkType::kSphere:
            return "sphere";
        default:
            return "unsupported type";
    }
}

}
}
