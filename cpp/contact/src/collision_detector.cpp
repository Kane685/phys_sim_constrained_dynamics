#include "contact/include/collision_detector.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "link/include/box_link.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

const CollisionEventGroup DetectPointCloudPlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation) {

    const std::string error_location = "contact::DetectPointCloudPlaneCollision";
    CheckCondition(first->type() == ContactProxyType::kPointCloud,
        error_location, "Incompatible first type.");
    CheckCondition(second->type() == ContactProxyType::kPlane,
        error_location, "Incompatible second type.");

    const auto point_cloud =
        std::dynamic_pointer_cast<PointCloudContactProxy>(first);
    const Isometry& point_cloud_transformation = first_transformation;
    const auto plane = std::dynamic_pointer_cast<PlaneContactProxy>(second);
    const Isometry& plane_transformation = second_transformation;

    CollisionEventGroup info(0);
    const Matrix3Xr world_point_cloud = ((point_cloud_transformation.first *
        point_cloud->vertices()).colwise() + point_cloud_transformation.second);
    // Convert the point to the plane frame.
    const Matrix3Xr plane_point_cloud = plane_transformation.first.transpose() *
        (world_point_cloud.colwise() - plane_transformation.second); 
    const RowVectorXr dists((plane->normal().transpose()
        * plane_point_cloud).array() + plane->offset());

    const Matrix3r world_plane_frame = plane_transformation.first *
        plane->local_frame();

    for (integer i = 0; i < point_cloud->vertex_num(); ++i) {
        if (dists(i) < 0) {
            // Create a contact event.
            CollisionEvent event;
            event.local_frame = world_plane_frame;
            event.collision_point = world_point_cloud.col(i);
            event.distance = dists(i);

            info.push_back(event);
        }
    }

    return info;
}

const CollisionEventGroup DetectSphereSphereCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation) {

    const std::string error_location = "contact::DetectSphereSphereCollision";
    CheckCondition(first->type() == ContactProxyType::kSphere,
        error_location, "Incompatible first type.");
    CheckCondition(second->type() == ContactProxyType::kSphere,
        error_location, "Incompatible second type.");

    const auto first_sphere =
        std::dynamic_pointer_cast<SphereContactProxy>(first);
    const Isometry& first_sphere_transformation = first_transformation;
    const auto second_sphere =
        std::dynamic_pointer_cast<SphereContactProxy>(second);
    const Isometry& second_sphere_transformation = second_transformation;

    CollisionEventGroup info(0);
    const Vector3r first_sphere_center = first_sphere_transformation.first *
        first_sphere->c() + first_sphere_transformation.second;
    const Vector3r second_sphere_center = second_sphere_transformation.first *
        second_sphere->c() + second_sphere_transformation.second;

    const real dist = (first_sphere_center - second_sphere_center).norm();
    const real r_sum = first_sphere->r() + second_sphere->r();
    if (dist < r_sum) {
        // Collision.
        CollisionEvent event;
        const Vector3r normal = first_sphere_center - second_sphere_center;
        const real normal_len = normal.stableNorm();
        CheckCondition(
            normal_len > std::max(first_sphere->r(), second_sphere->r()),
            error_location, "Spheres overlap too much.");
        event.local_frame = BuildFrameFromUnitNormal(
            normal.stableNormalized());
        event.collision_point = second_sphere_center +
            event.local_frame.col(2) *
            (second_sphere->r() - (r_sum - dist) / 2);
        event.distance = (dist - r_sum) / 2;

        info.push_back(event);
    }
    return info;
}

const CollisionEventGroup DetectSphereBoxCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation) {

    const std::string error_location = "contact::DetectSphereBoxCollision";
    CheckCondition(first->type() == ContactProxyType::kSphere,
        error_location, "Incompatible first type.");
    CheckCondition(second->type() == ContactProxyType::kBox,
        error_location, "Incompatible second type.");

    const auto sphere = std::dynamic_pointer_cast<SphereContactProxy>(first);
    const Isometry& sphere_transformation = first_transformation;
    const auto box = std::dynamic_pointer_cast<BoxContactProxy>(second);
    const Isometry& box_transformation = second_transformation;

    CollisionEventGroup info(0);

    ////////////////////////////////////////////////////////////////////////////
    // Task 3.2 (3 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Detect collision between a sphere and a box.
    //
    // BoxContactProxy creates an axis aligned box centered at the origin with a
    // size of size_. Next, it rotations the box with rotation_ followed by
    // translating it with translation_.
    //
    // Your implement consists of two parts:
    // - First, detect whether a collision event occurs. We suggest you
    //   convert the sphere's center to the box's body frame, which will make
    //   box-sphere collision detection easier.
    // - Next, if a collision event occurs, constructs the local frame, the
    //   collision point, and the penetration depth. The solution is not unique,
    //   and we encourage you to follow the lecture slides and/or the
    //   supplementary reading materials.
    //
    // TODO.

    // // center of sphere in world frame
    // Vector3r sphere_center = sphere_transformation.first*sphere->c()
    // +sphere_transformation.second;
    // // transform sphere center to box frame
    // Matrix3r box_world_R = box_transformation.first * box->R();
    // Vector3r box_world_t = box_transformation.second + box_transformation.first * box->t();
    // Matrix3r box_world_R_inv = box_world_R.transpose();
    // Vector3r sphere_center_box_local = box_world_R_inv * (sphere_center - box_world_t);
    // // collision detection
    // Vector3r half_size = box->s() * 0.5;
    // Vector3r closest_local;
    // closest_local.x() = std::max(-half_size.x(), std::min(sphere_center_box_local.x(), half_size.x()));
    // closest_local.y() = std::max(-half_size.y(), std::min(sphere_center_box_local.y(), half_size.y()));
    // closest_local.z() = std::max(-half_size.z(), std::min(sphere_center_box_local.z(), half_size.z()));
    // Vector3r closest_world = box_world_R * closest_local + box_world_t;
    // Vector3r diff = sphere_center - closest_world;
    // real dist_sq = diff.squaredNorm();
    // real sphere_radius = sphere->radius();
    // if (dist_sq < sphere_radius * sphere_radius) {
    //     real dist = std::sqrt(dist_sq);

    //     CollisionEvent event;
    //     // get collision point
    //     event.collision_point = closest_world;
    //     // get penetration distance
    //     event.distance = dist - sphere_radius;
    //     // get local frame
    //     Vector3r normal;
    //     normal = diff / dist;
    //     event.local_frame.col(2) = normal;
    //     event.local_frame.col(0) = event.local_frame.col(2).cross(Vector3r::UnitX()).normalized();
    //     event.local_frame.col(1) = event.local_frame.col(2).cross(event.local_frame.col(0));
    //     info.push_back(event);
    // }

    return info;
}

const CollisionEventGroup DetectSpherePlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation) {

    const std::string error_location = "contact::DetectSpherePlaneCollision";
    CheckCondition(first->type() == ContactProxyType::kSphere,
        error_location, "Incompatible first type.");
    CheckCondition(second->type() == ContactProxyType::kPlane,
        error_location, "Incompatible second type.");

    const auto sphere = std::dynamic_pointer_cast<SphereContactProxy>(first);
    const Isometry& sphere_transformation = first_transformation;
    const auto plane = std::dynamic_pointer_cast<PlaneContactProxy>(second);
    const Isometry& plane_transformation = second_transformation;

    CollisionEventGroup info(0);

    ////////////////////////////////////////////////////////////////////////////
    // Task 3.1 (3 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Detect collision between a sphere and an infinitely large plane and
    // return info, the collision information stored as a list of CollisionEvent
    // objects.
    //
    // The SphereContactProxy creates a sphere centered at center_ with a radius
    // radius_. sphere_transformation further moves this sphere to the world
    // frame with a rotation matrix sphere_transformation.first and a
    // translation vector sphere_transformation.second.
    //
    // Similarly, PlaneContactProxy creates a plane which is further transformed
    // by plane_transformation into the world frame. The plane represents a half
    // space with the normal pointing outward, i.e., the normal resides in the
    // exterior of the object represented by the infinitely large halfspace.
    //
    // Your implementation should determine whether the given sphere and plane
    // collides. If so, create a CollisionEvent object (one is enough) and
    // insert it into info.
    //
    // The CollisionEvent should create a local frame, record the point of
    // collision, and store the penetration depth.
    // - The local_frame is a 3 x 3 matrix whose three columns are orthonormal
    //   bases in the world frame. The z axis (the last column) should be the
    //   contact normal pointing from the plane to the sphere.
    // - The collision_point is a 3D point in the world frame. It is where we
    //   consider the collision occurs and the origin of the local_frame above.
    //   The choice is not unique. Here, we require that the point you return
    //   is on the surface of the sphere and that its offset with respect to the
    //   sphere's center is parallel to the contact normal in local_frame.
    // - Finally, the distance in CollisionEvent measures the penetration of the
    //   collision_point. We require that you implement the signed distance from
    //   collision_point to the plane and use a negative distance to indicate
    //   that the sphere collides with the plane.
    //
    // TODO.
    
    // center of sphere in world frame
    Vector3r sphere_center = sphere_transformation.first*sphere->c()
    +sphere_transformation.second;
    // normal of plane 
    Vector3r plane_normal = (second_transformation.first * plane->n()).normalized();
    Vector3r plane_point = plane_transformation.second;
    real signed_distance = (sphere_center - plane_point).dot(plane_normal);
    if (signed_distance < sphere->r()){
        CollisionEvent event;
        // get collision point
        event.collision_point = sphere_center - plane_normal*sphere->r();
        // geSt local frame
        event.local_frame.col(2) = (sphere_center-event.collision_point).normalized();
        event.local_frame.col(0) = event.local_frame.col(2).cross(Vector3r::UnitX()).normalized();
        event.local_frame.col(1) = event.local_frame.col(2).cross(event.local_frame.col(0));
        // get penetrative distance
        event.distance = signed_distance - sphere->r();
        info.push_back(event);
    }
        return info;
}

const CollisionEventGroup DetectBoxPlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation) {

    const std::string error_location = "contact::DetectBoxPlaneCollision";
    CheckCondition(first->type() == ContactProxyType::kBox,
        error_location, "Incompatible first type.");
    CheckCondition(second->type() == ContactProxyType::kPlane,
        error_location, "Incompatible second type.");

    const auto box = std::dynamic_pointer_cast<BoxContactProxy>(first);
    const Isometry& box_transformation = first_transformation;
    const auto plane = std::dynamic_pointer_cast<PlaneContactProxy>(second);
    const Isometry& plane_transformation = second_transformation;

    CollisionEventGroup info(0);
    // Collect vertices from the box.
    const Matrix3Xr unit_box_vertices = link::CreateBox().first;
    const integer vertex_num = static_cast<integer>(unit_box_vertices.cols());
    const Matrix3Xr box_vertices =
        (box->R() * unit_box_vertices.cwiseProduct(box->s() / 2 *
            RowVectorXr::Ones(vertex_num))).colwise() + box->t();
    const Matrix3Xr world_box_vertices =
        (box_transformation.first * box_vertices).colwise() +
        box_transformation.second;
    // Convert the point to the plane frame.
    const Matrix3Xr plane_box_vertices = plane_transformation.first.transpose()
        * (world_box_vertices.colwise() - plane_transformation.second);
    const RowVectorXr dists((plane->n().transpose()
        * plane_box_vertices).array() + plane->d());

    const Matrix3r world_plane_frame =
        plane_transformation.first * plane->local_frame();

    for (integer i = 0; i < vertex_num; ++i) {
        if (dists(i) < 0) {
            // Create a contact event.
            CollisionEvent event;
            event.local_frame = world_plane_frame;
            event.collision_point = world_box_vertices.col(i);
            event.distance = dists(i);

            info.push_back(event);
        }
    }

    return info;
}

void CollisionDetector::InitializeCollisionDetectionTable() {
    const integer cp_num = static_cast<integer>(ContactProxyType::kTotalNum);
    for (integer i = 0; i < cp_num; ++i)
        for (integer j = 0; j < cp_num; ++j)
            collision_detection_table_[i][j] = nullptr;

    // Register existing functions.
    const integer pc_idx = ToInt(ContactProxyType::kPointCloud);
    const integer sp_idx = ToInt(ContactProxyType::kSphere);
    const integer box_idx = ToInt(ContactProxyType::kBox);
    const integer pl_idx = ToInt(ContactProxyType::kPlane);
    collision_detection_table_[pc_idx][pl_idx] =
    collision_detection_table_[pl_idx][pc_idx] =
        DetectPointCloudPlaneCollisions;
    collision_detection_table_[sp_idx][sp_idx] =
        DetectSphereSphereCollisions;
    collision_detection_table_[sp_idx][box_idx] =
    collision_detection_table_[box_idx][sp_idx] =
        DetectSphereBoxCollisions;
    collision_detection_table_[sp_idx][pl_idx] =
    collision_detection_table_[pl_idx][sp_idx] =
        DetectSpherePlaneCollisions;
    collision_detection_table_[box_idx][pl_idx] =
    collision_detection_table_[pl_idx][box_idx] =
        DetectBoxPlaneCollisions;
}

const CollisionEventGroup CollisionDetector::DetectCollisions(
    const std::shared_ptr<ContactProxy>& first_contact_proxy,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second_contact_proxy,
    const Isometry& second_transformation) const {

    const std::string error_location =
        "contact::CollisionDetector::DetectCollision";

    CollisionEventGroup info;
    const integer cp_num = static_cast<integer>(ContactProxyType::kTotalNum);
    const integer first_idx = ToInt(first_contact_proxy->type());
    const integer second_idx = ToInt(second_contact_proxy->type());
    CheckCondition(0 <= first_idx && first_idx < cp_num &&
        0 <= second_idx && second_idx < cp_num, error_location,
        "Unsupported contact proxy type.");
    CheckCondition(collision_detection_table_[first_idx][second_idx],
        error_location, "Unsupported collision detection function.");
    if (first_idx <= second_idx) {
        info = collision_detection_table_[first_idx][second_idx](
            first_contact_proxy, first_transformation,
            second_contact_proxy, second_transformation);
    } else {
        // Swap order.
        info = collision_detection_table_[first_idx][second_idx](
            second_contact_proxy, second_transformation,
            first_contact_proxy, first_transformation);

        // Swap the two proxies.
        for (auto& event : info) {
            // Only the local frame needs to be swapped because it is assumed to
            // point from the second object to the first object.
            event.local_frame *= -1;
        }
    }
    return info;
}

}
}
