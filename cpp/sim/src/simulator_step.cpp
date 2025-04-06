#include "sim/include/simulator.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace sim {

void Simulator::ClearForce(const integer link_index) {
    CheckLinkIndex(link_index, "ClearForce");
    links_[link_index]->ClearForce();
}

void Simulator::ClearForces() {
    for (auto& link : links_) {
        link->ClearForce();
    }
}

void Simulator::ClearTorque(const integer link_index) {
    CheckLinkIndex(link_index, "ClearTorque");
    links_[link_index]->ClearTorque();
}

void Simulator::ClearTorques() {
    for (auto& link : links_) {
        link->ClearTorque();
    }
}

void Simulator::ApplyForce(const integer link_index,
    const Vector3r& attach_point, const Vector3r& force) {
    CheckLinkIndex(link_index, "ApplyForce");
    links_[link_index]->ApplyForce(attach_point, force);
}

void Simulator::ApplyTorque(const integer link_index, const Vector3r& torque) {
    CheckLinkIndex(link_index, "ApplyTorque");
    links_[link_index]->ApplyTorque(torque);
}

const StepInfo Simulator::Step(const real time_step, const Options& opt) {
    const std::string error_location = "sim::Simulator::Step";
    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real h = time_step;
    const real inv_h = ToReal(1) / h;

    // Compute inertia.
    std::vector<Matrix3r> inertia; inertia.reserve(link_num_);
    for (const auto& link : links_) inertia.push_back(link->ToWorldInertia());

    // Assemble the matrix.
    const SparseMatrixXr lhs = ComputeMomentumAndJointConstraintMatrix(
        h, inertia, opt);

    // Assemble the vector.
    const VectorXr rhs_m = ComputeMomentumVector(h, inertia, opt);
    const VectorXr rhs_j = ComputeJointConstraintVector(h, opt);
    VectorXr rhs = VectorXr::Zero(rhs_m.size() + rhs_j.size());
    rhs << rhs_m,
           rhs_j;

    // Collision detection.
    const ContactInfoGroup contacts = DetectCollisions();

    // Solve the system of equations.
    Matrix6Xr q_dot_plus;
    VectorXr mu_plus;
    if (contacts.empty()) {
        // Solve the system without any complementarity constraints.
        SolveDynamicsWithoutComplementarity(lhs, rhs, opt, q_dot_plus, mu_plus);
    } else {
        // Assemble the contact Jacobian.
        const SparseMatrixXr ct_J = ComputeContactJacobianMatrix(contacts);
        const VectorXr ct_neg_c_inv_h_sqr =
            ComputeContactComplianceVector(h, contacts);
        const VectorXr ct_neg_phi_inv_h =
            ComputeContactConstraintVector(h, contacts);
        SolveDynamicsWithComplementarity(lhs, rhs,
            contacts, ct_J, ct_neg_c_inv_h_sqr, ct_neg_phi_inv_h,
            opt, q_dot_plus, mu_plus);
    }

    // Write back results.
    StepLinksAndJoints(h, q_dot_plus, opt);

    // Collect info for debugging or visualization.
    StepInfo step_info;
    step_info.contacts = contacts;

    return step_info;
}

const SparseMatrixXr Simulator::ComputeMomentumAndJointConstraintMatrix(
    const real time_step, const std::vector<Matrix3r>& inertia,
    const Options& opt) const {

    const std::string error_location =
        "sim::Simulator::ComputeMomentumAndJointConstraintMatrix";

    const std::vector<Eigen::Triplet<real>> M_nonzeros = ToTriplet(
        ComputeMomentumMatrix(time_step, inertia, opt));
    const std::vector<Eigen::Triplet<real>> J_nonzeros = ToTriplet(
        ComputeJointConstraintJacobianMatrix(opt));
    const VectorXr D = ComputeJointConstraintComplianceVector(time_step, opt);

    // Assemble the matrix.
    std::vector<Eigen::Triplet<real>> nonzeros = M_nonzeros;
    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    for (const auto& ijv : J_nonzeros) {
        const integer row = 6 * link_num_ + ijv.row();
        const integer col = ijv.col();
        nonzeros.push_back(Eigen::Triplet<real>(row, col, ijv.value()));
        nonzeros.push_back(Eigen::Triplet<real>(col, row, ijv.value()));
    }
    for (integer i = 0; i < joints_dof; ++i)
        nonzeros.push_back(Eigen::Triplet<real>(6 * link_num_ + i,
            6 * link_num_ + i, D(i)));

    const SparseMatrixXr lhs = FromTriplet(6 * link_num_ + joints_dof,
        6 * link_num_ + joints_dof, nonzeros);

    return lhs;
}

const SparseMatrixXr Simulator::ComputeMomentumMatrix(const real time_step,
    const std::vector<Matrix3r>& inertia, const Options& opt) const {

    const std::string error_location =
        "sim::Simulator::ComputeMomentumMatrix";
    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real h = time_step;
    const real inv_h = ToReal(1) / h;

    // Assemble the matrix.
    std::vector<Eigen::Triplet<real>> nonzeros;

    ////////////////////////////////////////////////////////////////////////////
    // Task 2.2 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Assemble the mass matrix M used in the constrained dynamics. We write the
    // Newton-Euler equations from all links into the matrix form:
    //
    // Ma = f,
    //
    // where M is a (6 link_num_) x (6 link_num_) block-diagonal matrix, the
    // vector a includes all linear and angular accelerations from link_num_
    // links and is of dimension 6 link_num_, and f is a 6 link_num_-dimensional
    // vector. For Newton's equation, f is the net applied force; for Euler's
    // equation, we put everything other than I\dot{\omega} into f.
    //
    // This function represents M as a sparse matrix. Your job is to add all
    // nonzero entries in M into the std::vector nonzeros. You can access the
    // inertia tensor of the i-th link from the input inertia[i], which is
    // computed from Link::ToWorldInertia you implemented in Task 2.1. The
    // resultant M is a block-diagonal matrix. Each block is of size 6 x 6
    // corresponding to the mass matrix for one link's Newton-Euler equation.
    //
    // The sample code below shows how to insert a nonzero entry mij at the i-th
    // row and the j-th column of the matrix.
    //
    // TODO.
    // const integer i = 0;
    // const integer j = 0;
    // const real mij = 0;
    // nonzeros.emplace_back(i, j, mij);

    for (int i = 0; i < link_num_; ++i) {
        // start point of each M
        int start = i * 6;
    
        // Mass
        const real m = links_[i]->mass();
        // did not find mass in this file
        for (int d = 0; d < 3; ++d)
            nonzeros.emplace_back(start + d, start + d, m);
    
        // Inertia
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                nonzeros.emplace_back(start + 3 + r, start + 3 + c, inertia[i](r, c));
    }

    const SparseMatrixXr lhs = FromTriplet(6 * link_num_,
        6 * link_num_, nonzeros);

    return lhs;
}

const VectorXr Simulator::ComputeMomentumVector(const real time_step,
    const std::vector<Matrix3r>& inertia, const Options& opt) const {
    const std::string error_location =
        "sim::Simulator::ComputeMomentumVector";
    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real h = time_step;
    CheckCondition(static_cast<integer>(inertia.size()) == link_num_,
        error_location, "Wrong inertia size.");

    ////////////////////////////////////////////////////////////////////////////
    // Task 2.3 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Consider the temporal discretization of the matrix form Ma = f (See Task
    // 2.2):
    //
    // M (v+ - v) / h = f, or
    //
    // M v+ = hf + Mv.
    //
    // Here, v+ is the new velocity to be solved at the end of this time step
    // with a step size h. Your goal is to implement this function that returns
    // the right-hand side vector hf + Mv. The resultant vector contains
    // 6 link_num_ elements: the first 6 elements is the vector hf + Mv for the
    // first link, the next 6 elements for the second link, and so on. 
    //
    // TODO.

    VectorXr rhs = VectorXr::Zero(link_num_*6);
    for (int i=0;i<link_num_;++i){
        // define variables needed
        const real m = links_[i]->mass();
        const Vector3r& v = links_[i]->v();
        const Vector3r& w = links_[i]->omega();
        const Vector3r& f = links_[i]->f();
        const Vector3r& tau = links_[i]->tau();

        // calculate linear velocity and rotational velocity
        Vector3r linear_part = h*f + m*v;
        Vector3r angular_part = h*tau + inertia[i]*w;

        // get rhs
        rhs.segment<3>(i*6+0) = linear_part;
        rhs.segment<3>(i*6+3) = angular_part;
    }
    return rhs;
}

const SparseMatrixXr Simulator::ComputeJointConstraintJacobianMatrix(
    const Options& opt) const {

    const std::string error_location =
        "sim::Simulator::ComputeJointConstraintJacobianMatrix";

    // Assemble the matrix.
    std::vector<Eigen::Triplet<real>> nonzeros;
    const auto add_block = [&](const integer begin_row_idx,
        const integer begin_col_idx, const MatrixXr& matrix) -> void {
        const integer row_num = static_cast<integer>(matrix.rows());
        const integer col_num = static_cast<integer>(matrix.cols());
        for (integer i = 0; i < row_num; ++i)
            for (integer j = 0; j < col_num; ++j)
                nonzeros.push_back(Eigen::Triplet<real>(
                    begin_row_idx + i, begin_col_idx + j, matrix(i, j)));
    };
    for (const auto& joint : joints_) {
        // Add J.
        link::LinkGroup group;
        for (const integer i : joint.link_indices)
            group.push_back(links_[i]);
        const std::vector<MatrixX6r> Ji = joint.joint->Jphi(group);
        for (integer i = 0; i < joint.link_num; ++i) {
            const integer row_idx = joint.offset;
            const integer col_idx = 6 * joint.link_indices[i];
            add_block(row_idx, col_idx, Ji[i]);
        }
    }
    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    return FromTriplet(joints_dof, 6 * link_num_, nonzeros);
}

const VectorXr Simulator::ComputeJointConstraintComplianceVector(
    const real time_step, const Options& opt) const {

    const std::string error_location =
        "sim::Simulator::ComputeJointConstraintComplianceVector";
    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real h = time_step;
    const real inv_h = ToReal(1) / h;

    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    VectorXr compliance_vec = VectorXr::Zero(joints_dof);
    for (const auto& joint : joints_) {
        // Add -C / h2.
        const real val = -joint.joint->c() * inv_h * inv_h;
        compliance_vec.segment(joint.offset, joint.joint->m()) =
            VectorXr::Constant(joint.joint->m(), val);
    }

    return compliance_vec;
}

const VectorXr Simulator::ComputeJointConstraintVector(const real time_step,
    const Options& opt) const {
    const std::string error_location =
        "sim::Simulator::ComputeJointConstraintVector";
    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real inv_h = ToReal(1) / time_step;

    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    // Assemble -phi / h.
    VectorXr rhs = VectorXr::Zero(joints_dof);
    for (const auto& joint : joints_) {
        link::LinkGroup group;
        for (const integer i : joint.link_indices)
            group.push_back(links_[i]);
        // Note that we require joints to be equality constraints.
        rhs.segment(joint.offset, joint.joint->m())
            += -inv_h * joint.joint->phi(group);
    }
    return rhs;
}

void Simulator::SolveDynamicsWithoutComplementarity(
    const SparseMatrixXr& momentum_and_joint_constraint_matrix,
    const VectorXr& momentum_and_joint_constraint_vector,
    const Options& opt,
    Matrix6Xr& new_velocities, VectorXr& new_constraint_impulses) const {
    const std::string error_location =
        "sim::Simulator::SolveDynamicsWithoutComplementarity";

    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    const SparseMatrixXr& lhs = momentum_and_joint_constraint_matrix;
    CheckCondition(static_cast<integer>(lhs.rows()) == 6 * link_num_
        + joints_dof, error_location, "Incompatible lhs size.");
    const VectorXr& rhs = momentum_and_joint_constraint_vector;
    CheckCondition(rhs.size() == lhs.cols(), error_location,
        "Incompatible rhs size.");

    Eigen::SparseLU<SparseMatrixXr> sparse_solver;
    sparse_solver.compute(lhs);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Decomposition failed.");

    const VectorXr q_dot_and_mu_plus = sparse_solver.solve(rhs);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Solve failed.");

    new_velocities = q_dot_and_mu_plus.head(6 * link_num_).reshaped(
        6, link_num_);
    new_constraint_impulses = q_dot_and_mu_plus.tail(joints_dof);
}

void Simulator::SolveDynamicsWithComplementarity(
    const SparseMatrixXr& momentum_and_joint_constraint_matrix,
    const VectorXr& momentum_and_joint_constraint_vector,
    const ContactInfoGroup& contacts,
    const SparseMatrixXr& contact_jacobian_matrix,
    const VectorXr& contact_compliance_vector,
    const VectorXr& contact_constraint_vector,
    const Options& opt,
    Matrix6Xr& new_velocities, VectorXr& new_constraint_impulses) const {

    // The LCP problem:
    // [A,      J.T] [x ] = [b]
    // [J,      D  ] [mu]   [v].
    // =>
    // [J,   JA^-1 J.T] [x ] = [JA^-1 b]
    // [J,   D        ] [mu] = [v].
    // (JA^-1 J.T - D) mu = JA^-1 b - v.
    // =>
    // v = (JA^-1 J.T - D)(-mu) + JA^-1 b.

    const std::string error_location =
        "sim::Simulator::SolveDynamicsWithComplementarity";

    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }
    const integer contact_num = static_cast<integer>(contacts.size());
    CheckCondition(contact_num > 0, error_location,
        "Empty commplemnetarity constraints.");

    const SparseMatrixXr& A = momentum_and_joint_constraint_matrix;
    CheckCondition(A.rows() == A.cols() &&
        static_cast<integer>(A.rows()) == 6 * link_num_ + joints_dof,
        error_location, "Incorrect A size.");
    const VectorXr& b = momentum_and_joint_constraint_vector;
    CheckCondition(b.size() == A.cols(), error_location, "Incorrect b size.");
    const SparseMatrixXr& Jc = contact_jacobian_matrix;
    CheckCondition(Jc.cols() == link_num_ * 6 &&
        static_cast<integer>(Jc.rows()) == contact_num * 3, error_location,
        "Incorrect Jc size.");
    const std::vector<Eigen::Triplet<real>> Jc_nonzeros = ToTriplet(Jc);
    std::vector<Eigen::Triplet<real>> Jt_nonzeros;
    Jt_nonzeros.reserve(Jc_nonzeros.size());
    for (const auto& triplet : Jc_nonzeros)
        Jt_nonzeros.push_back(Eigen::Triplet<real>(
            triplet.col(), triplet.row(), triplet.value()));
    const SparseMatrixXr Jt = FromTriplet(6 * link_num_ + joints_dof,
        contact_num * 3, Jt_nonzeros);

    CheckCondition(static_cast<integer>(contact_compliance_vector.size()) ==
        contact_num * 3, error_location, "Incorrect compliance vector size.");
    const VectorXr& D = contact_compliance_vector;

    // v = M * (-mu) + N.
    // M = (J A^-1 J.T - D).
    // N = J A^-1 b.

    // Factorize A.
    Eigen::SparseLU<SparseMatrixXr> sparse_solver;
    sparse_solver.compute(A);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Decomposition failed.");

    // Solve A^-1 b.
    const VectorXr A_inv_b = sparse_solver.solve(b);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Solve b failed.");
    const VectorXr N(A_inv_b.transpose() * Jt);

    // Solve A^-1 J.T.
    const MatrixXr A_inv_Jt = sparse_solver.solve(Jt);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Solve Jt failed.");
    MatrixXr M = Jt.transpose() * A_inv_Jt;
    // Subtract diagonal D.
    for (integer i = 0; i < contact_num * 3; ++i)
        M(i, i) -= D(i);

    // The block-LCP algorithm:
    // v = M * (-mu) + N.
    const integer max_blcp_iter = opt.integer_option().at("max_blcp_iter");
    CheckCondition(max_blcp_iter > 0, error_location, "Invalid max_blcp_iter.");

    // For contact:
    // Our construction uses neg_mu to represent the contact and friction
    // impulse given by the second object to the first object in the contact.
    // The impulse is represented in the local frame whose normal points from
    // the second object to the first one.
    // - If both objects are links, then J.T * neg_mu represents the impulse
    //   each of the object experiences.
    // - If only the first object is a link, then J.T * neg_mu is the impulse
    //   it experiences from the environmental obstacle.
    // - If only the second object is a link, then J.T * neg_mu is the impulse
    //   it experiences from the environmental obstacle (the first object in our
    //   collision detection in this case).
    // In summary, regardless of the nature of the contact (link vs. link or
    // link vs env.), we can consistently interpret neg_mu as the impulse from
    // the second object to the first object and require its normal component
    // to be nonnegative. The same is true for relative velocities. There is no
    // need to flip signs based on the order of link vs. env. or whatsoever.
    VectorXr neg_mu = VectorXr::Zero(contact_num * 3);
    std::vector<Matrix3r> Mii(0); Mii.reserve(contact_num);
    std::vector<Eigen::HouseholderQR<Matrix3r>>Mii_inverse(0);
    Mii_inverse.reserve(contact_num);
    for (integer i = 0; i < contact_num; ++i) {
        Mii.push_back(M.block<3, 3>(3 * i, 3 * i));
        const Eigen::HouseholderQR<Matrix3r> qr(Mii.back());
        Mii_inverse.push_back(qr);
    }
    VectorXr last_neg_mu = neg_mu;
    bool converged = false;
    const real impulse_tol = opt.real_option().at("blcp_impulse_abs_tol");
    CheckCondition(impulse_tol > 0, error_location, "The impulse tolerance "
        "should be strictly positive.");
    for (integer k = 0; k < max_blcp_iter; ++k) {
        last_neg_mu = neg_mu;
        // Contacts.
        for (integer i = 0; i < contact_num; ++i) {
            const Vector3r neg_mui = neg_mu.segment<3>(i * 3);
            const Vector3r Ni = N.segment<3>(i * 3);
            // Freeze other blocks and solve the 3 x 3 LCP.
            // v[i] = M[ii] * neg_mu[i] + M[ij] * neg_mu[j] + N[i].
            // Here, neg_mu is the contact and friction impulse in the local
            // frame. Taking the 3D case as a concrete example, neg_mu.x and
            // neg_mu.y are the frictions and neg_mu.z is the contact impulse.
            //
            // The resultant v[i] refers to the relative velocity of the first
            // object w.r.t. the second object. The compliance stabilization is
            // applied to the normal velocity only.
            //
            // The BLCP condition:
            // - Contact:
            //   * v.z = -phi / h, neg_mu.z > 0.
            //   * v.x > 0 => neg_mu.x = -mu * neg_mu.z
            //     v.x < 0 => neg_mu.x = mu * neg_mu.z.
            //     v.x = 0 => -mu * neg_mu.z <= neg_mu.x <= mu * neg_mu.z.
            // - No contact:
            //   * neg_mu == 0 && v.z >= -phi / h.
            const Vector3r b = M.middleRows<3>(i * 3) * neg_mu
                - Mii[i] * neg_mui + Ni;
            // v[i] = M[ii] * neg_mu[i] + b.
            const real neg_phi_inv_h = contact_constraint_vector(i * 3 + 2);
            neg_mu.segment<3>(i * 3) = SolveLocalFrictionalBlcp(Mii[i],
                Mii_inverse[i], b, neg_phi_inv_h,
                contacts[i].frictional_coefficient);
        }
        // Convergence check.
        if ((neg_mu - last_neg_mu).cwiseAbs().maxCoeff() <= impulse_tol) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        PrintWarning(error_location, "BLCP failed to converge.");
    }

    // Write back results.
    // A * x + J.T * mu = b.
    const VectorXr x = sparse_solver.solve(b + Jt * neg_mu);
    CheckCondition(sparse_solver.info() == Eigen::Success,
        error_location, "Solve b + Jt * neg_mu failed.");
    new_velocities = x.head(6 * link_num_).reshaped(6, link_num_);
    new_constraint_impulses = x.tail(joints_dof);
}

void Simulator::StepLinksAndJoints(const real time_step,
    const Matrix6Xr& new_velocities, const Options& opt) {

    const std::string error_location =
        "sim::Simulator::StepLinksAndJoints";

    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real h = time_step;
    CheckCondition(static_cast<integer>(new_velocities.cols()) == link_num_,
        error_location, "Incompatible new_velocities size.");

    integer joints_dof = 0;
    if (!joints_.empty()) {
        joints_dof = joints_.back().offset + joints_.back().joint->m();
    }

    // Write back results.
    for (integer i = 0; i < link_num_; ++i) {
        const auto& link = links_[i];
        link->set_v(new_velocities.col(i).head(3));
        link->set_omega(new_velocities.col(i).tail(3));
        link->set_t(link->t() + h * link->v());
        link->set_R(ToRotationMatrix(link->omega() * h) * link->R());
    }
}

}
}
