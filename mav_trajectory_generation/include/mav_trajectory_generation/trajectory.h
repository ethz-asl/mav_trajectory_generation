#ifndef MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_
#define MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_

namespace mav_trajectory_generation {

class Trajectory {
 public:
  Trajectory();

  Eigen::VectorXd evaluate(double t, int derivative_order) const;

  int D() const { return D_; }
  int N() const { return N_; }
  int K() const { return segments_.size(); }


 private:
  int D_;  // Number of dimensions.
  int N_;  // Number of coefficients.

  // TODO(helenol): align!!!!!!!!11!1!
  // K is number of segments...
  std::vector<Segment> segments_;
};

}  // namespace mav_trajectory_generation


#endif  // MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_
