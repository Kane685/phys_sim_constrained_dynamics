#include "joint/include/binary_hinge_joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

const VectorXr BinaryHingeJoint::phi(const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "BinaryHingeJoint::phi");

    ////////////////////////////////////////////////////////////////////////////
    // Task 1.2 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Finish the definition of the constraint function for a hinge joint
    // connecting two links, q[0] and q[1].
    //
    // The joint attaches first_attach_point_ in q[0] to second_attach_point_
    // in q[1]. Both first_attach_point_ and second_attach_point_ are 3D
    // coordinates of the points in the body frames of q[0] and q[1],
    // respectively.
    //
    // q[0] and q[1] spin around a shared axis passing first_attach_point_ (
    // second_attach_point_). The axis' direction is first_attach_direction_
    // when expressed as a vector in the body frame of q[0] and
    // second_attach_direction_ in the body frame of q[1].
    //
    // Your goal is to implement a 6-dimensional constraint function
    // phi(q[0], q[1]) so that phi = 0 imposes the above conditions.
    //
    // Reading the phi functions in other joint classes may be helpful.
    //
    // TODO.

    // set empty 6-dimension vector
    Vector6r phi; phi.setZero();
    // define constraint function
    phi.head(3) = q[0]->ToWorldPoint(first_attach_point_)
                  - q[1]->ToWorldPoint(second_attach_point_);

    phi.tail(3) = q[0]->ToWorldPoint(first_attach_point_+first_attach_direction_)
                  - q[1]->ToWorldPoint(second_attach_point_+second_attach_direction_);


    return phi;
}


const std::vector<MatrixX6r> BinaryHingeJoint::Jphi(
    const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "BinaryHingeJoint::Jphi");

    std::vector<MatrixX6r> J(this->n(), MatrixXr::Zero(this->m(), 6));

    ////////////////////////////////////////////////////////////////////////////
    // Task 1.3 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Implement the Jacobian of the phi function above so that
    //
    // -J[i].T * compliance_^-1 * phi (i = 0 or 1) describes the constraint
    // force induced from the quadratic energy V = compliance_ * phi^2 / 2.
    // Check the lecture slides for more details. Reading the Jphi functions in
    // other joint classes can also be helpful.
    //
    // TODO.

    J[0].topRows(3) = q[0]->ComputePointJacobian(first_attach_point_);
    J[0].bottomRows(3) = q[0]->ComputeVectorJacobian(first_attach_point_+first_attach_direction_);
    J[1].topRows(3) = -q[1]->ComputePointJacobian(second_attach_point_);
    J[1].bottomRows(3) = -q[1]->ComputeVectorJacobian(second_attach_point_+second_attach_direction_);
    
    return J;
}

void BinaryHingeJoint::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "joint::BinaryHingeJoint::InitializeDerived";

    for (const std::string& key : { "first_attach_point", "second_attach_point",
        "first_attach_direction", "second_attach_direction" }) {
        CheckCondition(opt.GetVectorOptionSize(key) == 3, error_location,
            "Incompatible " + key + " size.");
    }
    first_attach_point_ = opt.vector_option().at("first_attach_point");
    second_attach_point_ = opt.vector_option().at("second_attach_point");
    first_attach_direction_ =
        opt.vector_option().at("first_attach_direction").stableNormalized();
    second_attach_direction_ =
        opt.vector_option().at("second_attach_direction").stableNormalized();
    CheckCondition(IsClose(first_attach_direction_.stableNorm(), 1,
        ToReal(1e-6), ToReal(1e-3)), error_location,
        "Too small first_attach_direction.");
    CheckCondition(IsClose(second_attach_direction_.stableNorm(), 1,
        ToReal(1e-6), ToReal(1e-3)), error_location,
        "Too small second_attach_direction.");
}

}
}
