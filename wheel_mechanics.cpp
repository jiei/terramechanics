#include "wheel_mechanics.h"

void Wheel_Mechanics::InputParam(soil input_soil, wheel input_wheel) {
  regolith = input_soil;
  rover_wheel = input_wheel;
  //パラメータ表示
  std::cout << "Parameters" << '\n';
  std::cout << "\tregolith" << '\n';
  std::cout << "\t\tc\t=\t" << regolith.c << '\n';
  std::cout << "\t\tk_c\t=\t" << regolith.k_c << '\n';
  std::cout << "\t\tphi\t=\t" << regolith.phi << '\n';
  std::cout << "\t\tk_phi\t=\t" << regolith.k_phi << '\n';
  std::cout << "\t\tn\t=\t" << regolith.n << '\n';
  std::cout << "\t\ta_0\t=\t" << regolith.a_0 << '\n';
  std::cout << "\t\ta_1\t=\t" << regolith.a_1 << '\n';
  std::cout << "\t\tk_x\t=\t" << regolith.k_x << '\n';
  std::cout << "\t\tPHI\t=\t" << regolith.PHI << '\n';
  std::cout << "\twheel" << '\n';
  std::cout << "\t\tW\t=\t" << rover_wheel.W << '\n';
  std::cout << "\t\tr\t=\t" << rover_wheel.r << '\n';
  std::cout << "\t\tb\t=\t" << rover_wheel.b << '\n';
  std::cout << "\t\tm\t=\t" << rover_wheel.m << '\n';
  std::cout << "\t\tI\t=\t" << rover_wheel.I << '\n';
  std::cout << "\t\tkai\t=\t" << rover_wheel.kai << '\n';
}

void Wheel_Mechanics::CalcStaticSinkage() {
  double a = 0, b = M_PI / 2, s;
  double eps = 0.0001;
  double tmp = pow(rover_wheel.r, regolith.n + 1) *
               (regolith.k_c + regolith.k_phi * rover_wheel.b) /
               pow(cos(regolith.PHI), regolith.n);

  std::cout << "Calculating Static Sinkage.........." << '\n';

  while (fabs(a - b) > eps) {
    s = (a + b) / 2.0;
    if ((rover_wheel.W - tmp * StaticSinkageIntegration(s)) *
            (rover_wheel.W - tmp * StaticSinkageIntegration(a)) <
        0)
      b = s;
    else
      a = s;
  }
  h_s = rover_wheel.r * (1 - cos((a + b) / 2 - regolith.PHI));

  //#ifdef DEBUG
  std::cout << "----------Calculated Static Sinkage----------" << '\n';
  std::cout << "a = " << a << '\n';
  std::cout << "b = " << b << '\n';
  std::cout << "h_s = " << h_s << '\n';
  //#endif
}

inline double Wheel_Mechanics::StaticSinkageIntegration(double theta_s) {
  double theta = (-1) * theta_s + 2 * regolith.PHI, S = 0;

  while (theta < theta_s) {
    S += (StaticSinkageIntegratedFunc(theta, theta_s) +
          StaticSinkageIntegratedFunc(theta + d_theta, theta_s)) *
         d_theta / 2.0;
    theta += d_theta;
  }

  return S;
}

void Wheel_Mechanics::CalcWheelDynamics() {
  double time_tmp = 0.0;
  double x_2dot_tmp = 0.0, x_dot_tmp = 0.0, x_tmp = 0.0;
  double z_2dot_tmp = 0.0, z_dot_tmp = 0.0, z_tmp = 0.0;
  double h_tmp = h_s;

  while (time_tmp < time_limit) {
    time.push_back(time_tmp);

    x_2dot.push_back(x_2dot_tmp), x_dot.push_back(x_dot_tmp),
        x.push_back(x_tmp);
    z_2dot.push_back(z_2dot_tmp), z_dot.push_back(z_dot_tmp),
        z.push_back(z_tmp);
    h.push_back(h_tmp);

    slip_ratio.push_back(CalcSlipRatio(x_dot.back()));

    theta_f = CalcTheta_f(h.back());
    theta_r = CalcTheta_r(h.back());
    theta_m = CalcTheta_m(slip_ratio.back());

    F_Integration(&F_x, &F_R, &F_z);
    F_DP.push_back(F_x.back() - F_R.back());

    x_2dot_tmp =
        (4 * F_DP.back() - ROVER_MASS * G * sin(regolith.PHI)) / ROVER_MASS;
    // x_2dot_tmp = F_DP.back() / rover_wheel.m;
    x_dot_tmp += x_2dot_tmp * dT;
    x_tmp += x_dot_tmp * dT;

    z_2dot_tmp = G * cos(regolith.PHI) - 4 * F_z.back() / ROVER_MASS;
    // z_2dot_tmp = G - F_z.back() / rover_wheel.m;
    z_dot_tmp += z_2dot_tmp * dT;
    z_tmp += z_dot_tmp * dT;

    h_tmp = h_s + z_tmp;

    time_tmp += dT;
  }

  std::cout << "----------Calculated Wheel Dynamics----------" << '\n';
  std::cout << "Status at last place" << '\n';
  std::cout << "time\t=\t" << time.back() << '\n';
  std::cout << "r_alpha_dot\t=\t" << r_alpha_dot << '\n';
  std::cout << "x\t=\t" << x.back() << "\tsize\t=\t" << x.size() << '\n';
  std::cout << "z\t=\t" << z.back() << "\tsize\t=\t" << z.size() << '\n';
  std::cout << "x_dot\t=\t" << x_dot.back() << "\tsize\t=\t" << x_dot.size()
            << '\n';
  std::cout << "z_dot\t=\t" << z_dot.back() << "\tsize\t=\t" << z_dot.size()
            << '\n';
  std::cout << "h\t=\t" << h.back() << "\tsize\t=\t" << h.size() << '\n';
  std::cout << "slip_ratio\t=\t" << slip_ratio.back() << "\tsize\t=\t"
            << slip_ratio.size() << '\n';
}

void Wheel_Mechanics::InitializeState(double dT_, double time_limit_) {
  dT = dT_;
  time_limit = time_limit_;
}

inline double Wheel_Mechanics::CalcSlipRatio(double x_dot_) {
  double slip_ratio_tmp =
      (r_alpha_dot - x_dot_) / r_alpha_dot; //ドライブ状態のみ考慮

  if (slip_ratio_tmp < -1) {
    std::cout << "slip ratio is below -1.0 !" << '\n';
    return -1.0;
  } else if (slip_ratio_tmp > 1) {
    std::cout << "slip ratio is over 1.0 !" << '\n';
    return 1.0;
  } else
    return slip_ratio_tmp;
}

inline double Wheel_Mechanics::CalcTheta_f(double h_) {
  return acos(1 - h_ / rover_wheel.r);
}

inline double Wheel_Mechanics::CalcTheta_r(double h_) {
  return (-1) * acos(1 - h_ * rover_wheel.kai / rover_wheel.r);
}

inline double Wheel_Mechanics::CalcTheta_m(double slip_ratio_) {
  return (regolith.a_0 + regolith.a_1 * slip_ratio_) * theta_f;
}

inline double Wheel_Mechanics::CalcJ_theta(double theta, double slip_ratio_) {
  return rover_wheel.r *
         (theta_f - theta - (1 - slip_ratio_) * (sin(theta_f) - sin(theta)));
}

inline double Wheel_Mechanics::CalcSigma_theta(double theta) {
  if (theta_m <= theta && theta <= theta_f) {
    return (regolith.k_c / rover_wheel.b + regolith.k_phi) *
               pow(rover_wheel.r * (cos(theta) - cos(theta_f)), regolith.n) +
           30000 * z_dot.back();
  } else if (theta_r <= theta && theta <= theta_m) {
    return (regolith.k_c / rover_wheel.b + regolith.k_phi) *
               pow(rover_wheel.r *
                       (cos(theta_f - (theta - theta_r) * (theta_f - theta_m) /
                                          (theta_m - theta_r)) -
                        cos(theta_f)),
                   regolith.n) +
           30000 * z_dot.back();
  } else {
    std::cout << "theta does not match the range @CalcSigma_theta" << '\n';
    return 0;
  }
}

inline double Wheel_Mechanics::CalcTau_theta() {
  return (regolith.c + sigma_theta * tan(regolith.phi)) *
         (1 - exp((-1) * j_theta / regolith.k_x));
}

inline bool Wheel_Mechanics::F_Integration(std::vector<double> *F_x_,
                                           std::vector<double> *F_R_,
                                           std::vector<double> *F_z_) {
  double theta = theta_r;
  double S_x = 0, S_R = 0, S_z = 0;

  while (theta < theta_f) {
    j_theta = CalcJ_theta(theta, slip_ratio.back());
    sigma_theta = CalcSigma_theta(theta);
    tau_theta = CalcTau_theta();

    S_x += F_xIntegratedFunc(theta) * d_theta;
    S_R += F_RIntegratedFunc(theta) * d_theta;
    S_z += F_zIntegratedFunc(theta) * d_theta;
    theta += d_theta;
  }

  F_x_->push_back(rover_wheel.r * rover_wheel.b * S_x);
  F_R_->push_back(rover_wheel.r * rover_wheel.b * S_R);
  F_z_->push_back(rover_wheel.r * rover_wheel.b * S_z);

  return true;
}
