/**
 *
 * @file main.cpp
 * @brief Main file for numerical simulation of wheel behavior on regolith
 * surface.
 * @author Takanishi Labratory, Environmental Monitoring team, Jiei Suzuki.
 * @date 2018/7/14
 *
 */

#define _USE_MATH_DEFINES
//#define SAVE_AS_CSV

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "matplotlib-cpp-master/matplotlibcpp.h" //グラフの描画
#include "wheel_mechanics.h"                     //メインのクラス

using namespace std;
namespace plt = matplotlibcpp;

/**
 * プログラムの簡単な流れ
 * 1.土壌・車輪パラメータの入力
 * 2.静的沈下量の算出
 * 3.動的変化量の逐次計算（オイラー法）
 * 4.グラフへの動的変化量の描画
 * 5.CSVファイルへの出力
 */

int main(int argc, char const *argv[]) {
/* code */
#ifdef SAVE_AS_CSV
  if (argc < 2) {
    std::cerr << "Invalid argment : Must input output-file name" << '\n';
  }
  std::string fname = "/home/monitoring/source/terramechanics/OutputData/" +
                      std::string(argv[1]);
  std::ofstream fout(fname.c_str());
  if (!fout) {
    std::cout << "Cannot open " << fname.c_str() << '\n';
    return -1;
  } else {
    std::cout << "Opened " << fname.c_str() << '\n';
  }
#endif

  Wheel_Mechanics w(0.1);

  soil regolith = {0.8,  1370,  M_PI * 37.2 / 180, 814000, 1.0, 0.4,
                   0.15, 0.025, M_PI * 10 / 180};
  wheel rover_wheel = {ROVER_MASS * G / 4, 0.15, 0.1, 6.5, 0.065, 1.1};

  w.InputParam(regolith, rover_wheel);
  w.CalcStaticSinkage();
  w.InitializeState(0.001, 3.0);
  w.CalcWheelDynamics();

  //簡易的にグラフで表示
  plt::plot(w.time, w.slip_ratio);
  plt::title("slip ratio");
  plt::xlabel("time t[s]");
  plt::ylabel("slip_ratio s");
  plt::show();
  plt::plot(w.time, w.h);
  plt::title("sinkage");
  plt::xlabel("time t[s]");
  plt::ylabel("sinkage h[m]");
  plt::show();
  plt::plot(w.time, w.x_dot);
  plt::title("forward velosity");
  plt::xlabel("time t[s]");
  plt::ylabel("forward_velosity dx/dt[m/s]");
  plt::show();
  plt::plot(w.time, w.F_DP);
  plt::title("traction force");
  plt::xlabel("time t[s]");
  plt::ylabel("traction force F_DP[N]");
  plt::show();

#ifdef SAVE_AS_CSV
  // CSVファイルに出力
  std::cout << "----------Saving to " << argv[1] << "----------" << '\n';
  fout << "time" << ',' << "x_2dot" << ',' << "x_dot" << ',' << "x" << ','
       << "z_2dot" << ',' << "z_dot" << ',' << "z" << ',' << "slip_ratio" << ','
       << "F_DP" << ',' << "F_z" << ',' << "h" << std::endl;
  for (int i = 0; i < w.time.size(); i++) {
    fout << w.time.at(i) << ',' << w.x_2dot.at(i) << ',' << w.x_dot.at(i) << ','
         << w.x.at(i) << ',' << w.z_2dot.at(i) << ',' << w.z_dot.at(i) << ','
         << w.z.at(i) << ',' << w.slip_ratio.at(i) << ',' << w.F_DP.at(i) << ','
         << w.F_z.at(i) << ',' << w.h.at(i) << std::endl;
  }
  fout.close();
#endif

  return 0;
}
