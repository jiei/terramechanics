#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

#define DEBUG

#define G 1.622       //月の重力
#define ROVER_MASS 10 //ローバ重量

//パラメータ
typedef struct {
  double c;     //土壌粘着力
  double k_c;   // cに依存する土壌変形定数
  double phi;   //内部摩擦角
  double k_phi; // phiに依存する土壌変形定数
  double n;     //指数定数
  double a_0;
  double a_1;
  double k_x; //車輪-土壌パラメータ
  double PHI; //斜面角度
} soil;       //土壌定数

typedef struct {
  double W;   //車輪にかかる荷重
  double r;   //車輪半径
  double b;   //車輪幅
  double m;   //車輪質量
  double I;   //慣性モーメント
  double kai; //離脱角係数
} wheel;      //車輪定数

class Wheel_Mechanics {
public:
  //土壌と車輪のパラメータ
  soil regolith;
  wheel rover_wheel;

  double dT;
  double time_limit;
  std::vector<double> time;

  double d_theta;

  //状態変数
  double r_alpha_dot;
  std::vector<double> x_2dot, x_dot, x;
  std::vector<double> z_2dot, z_dot, z;

  std::vector<double> slip_ratio;
  std::vector<double> F_x;
  std::vector<double> F_R;
  std::vector<double> F_DP;
  std::vector<double> F_z;
  std::vector<double> h;

  double h_s; //静的沈下量
  double sigma_theta;
  double tau_theta;
  double theta_m, theta_f, theta_r;
  double j_theta;

  Wheel_Mechanics(double r_alpha_dot_) { //コンストラクタ 引数は車輪速度
    r_alpha_dot = r_alpha_dot_;
    d_theta = 0.0001; //計算刻み幅
  }

  void InputParam(soil input_soil,
                  wheel input_wheel); //土壌・車輪パラメータ入力

  //静的沈下量を計算するのに必要な関数群
  void CalcStaticSinkage(); //静的沈下量計算
  inline double StaticSinkageIntegration(double theta_s);
  inline double StaticSinkageIntegratedFunc(double theta, double theta_s) {
    return pow((cos(theta - regolith.PHI) - cos(theta_s - regolith.PHI)),
               regolith.n) *
           cos(theta);
  }

  void CalcWheelDynamics(); //車輪ダイナミクス計算のメイン関数
  void InitializeState(double dT_, double time_limit_); //もろもろの変数を初期化
  inline double CalcSlipRatio(double x_dot_);    //スリップ率計算
  inline double CalcTheta_f(double h_);          //θfの計算
  inline double CalcTheta_r(double h_);          //θrの計算
  inline double CalcTheta_m(double slip_ratio_); //θmの計算
  inline double CalcJ_theta(double theta, double slip_ratio_); //ｊ（θ）の計算
  inline double CalcSigma_theta(double theta); //σ（θ）の計算
  inline double CalcTau_theta();               //τ(θ)の計算

  inline bool
  F_Integration(std::vector<double> *F_x_, std::vector<double> *F_R_,
                std::vector<double> *F_z_); //θで積分する操作をまとめた

  inline double F_xIntegratedFunc(double theta) {
    return tau_theta * cos(theta);
  }
  inline double F_RIntegratedFunc(double theta) {
    return sigma_theta * sin(theta);
  }
  inline double F_zIntegratedFunc(double theta) {
    return sigma_theta * cos(theta) + tau_theta * sin(theta);
  }
};
