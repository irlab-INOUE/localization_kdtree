#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <opencv2/opencv.hpp>

#include "KDTree.h"
#include "Map.h"

#define DEBUG_DISPLAY

struct UrgData {
  std::string type;
  long long timestamp1;
  long long timestamp2;
  int size;
  double start_angle;
  double end_angle;
  double step;
  int echo_size;
  std::vector<long> dat;
};

class Match : public Point {
public:
  double heading;
  double match;

  Match() {};
  Match(double x, double y, double a, double val) {
    this->x = x;
    this->y = y;
    this->heading = a;
    this->match = val;
  }
};

double evaluate(double rx, double ry, double ra, UrgData &data, int skip, MapClass &map, Node *root, cv::Mat &img) {
  double cs = cos(ra);
  double sn = sin(ra);
#ifdef DEBUG_DISPLAY
  cv::Mat display;
  img.copyTo(display);
#endif
  double dist = 0;
  for (int i = 0; i < data.dat.size(); i+=skip) {
    if (data.dat[i] > 35000) continue;
    double x = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
    double y = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

    // 座標変換
    double cx = cs * x - sn * y + rx;
    double cy = sn * x + cs * y + ry;

    // 探索する点
    Point addPoint(cx, cy);
    int ax =  addPoint.x / map.csize + map.ORIGIN_X;
    int ay =-(addPoint.y / map.csize) + map.ORIGIN_Y;
#ifdef DEBUG_DISPLAY
    cv::circle(display, cv::Point(ax, ay), 2, cv::Scalar(0, 200, 0), -1, cv::LINE_8);
#endif

    // 最近傍点探索
    Node *nns = nearest_neighbor_search_KDTree(addPoint, root, 0);
    int nx =  nns->pt.x / map.csize + map.ORIGIN_X;
    int ny =-(nns->pt.y / map.csize) + map.ORIGIN_Y;
#ifdef DEBUG_DISPLAY
    cv::circle(display, cv::Point(nx, ny), 4, cv::Scalar(200, 0, 0), -1, cv::LINE_8);
#endif

    // 点間距離
    dist += std::hypot(cx - nns->pt.x, cy - nns->pt.y);
    // 対応点
#ifdef DEBUG_DISPLAY
    cv::line(display, cv::Point(ax, ay), cv::Point(nx, ny), cv::Scalar(100, 100, 250), 1, cv::LINE_8);
#endif
  }
  dist /= static_cast<double>(data.dat.size());
  // ロボット
#ifdef DEBUG_DISPLAY
  cv::circle(display, cv::Point(rx/map.csize + map.ORIGIN_X, -ry/map.csize + map.ORIGIN_Y), 
             0.5/map.csize, cv::Scalar(200, 0, 200), 2);
  cv::imshow("RESULT", display);
  cv::waitKey(1);
#endif
  return dist;
}

int main(int argc, char *argv[]) {
  double xmin = -5.0;
  double xmax =  5.0;
  double ymin = -5.0;
  double ymax =  5.0;
  double dd = 0.1;
  int skip = 20;        // lidarのスキップ数

  std::vector<Match> match;
  MapClass map;
  cv::Mat img_map = cv::imread("../occMap.png");
  map.img2map(img_map, 0.025, 4, 82+50, 322+50);

  std::ifstream fin("../urglog_oneshot");

  UrgData data;
  fin >> data.type
    >> data.timestamp1
    >> data.size
    >> data.start_angle
    >> data.end_angle
    >> data.step
    >> data.echo_size;

  for (int i = 0; i < data.size/data.echo_size; i++) {
    long tmp;
    fin >> tmp;
    data.dat.emplace_back(tmp);
    fin >> tmp >> tmp;
  }
  fin >> data.timestamp2;

  std::vector<Point> XYs;
  for (int y = 0; y < map.HEIGHT; y++) {
    for (int x = 0; x < map.WIDTH; x++) {
      int check = map.get_occ(x, y);
      if (check < 20) {
        XYs.emplace_back((x - map.ORIGIN_X) * map.csize, 
                         (-y + map.ORIGIN_Y) * map.csize);
      }
    }
  }
  // KD-Tree の作成
  int depth = 0; 
  int start_index = 0;
  int last_index = XYs.size();
  Node *root_parent_is_nullptr = nullptr;
  Node *root = makeKDTree(XYs, start_index, last_index, depth, root_parent_is_nullptr);

  // 結果の表示
  cv::Mat img = cv::Mat(cv::Size(map.WIDTH, map.HEIGHT), CV_8UC3, cv::Scalar(182, 182, 182));
  cv::line(img, cv::Point(0, map.ORIGIN_Y), cv::Point(map.WIDTH, map.ORIGIN_Y), cv::Scalar(80, 80, 80), 1, cv::LINE_8);
  cv::line(img, cv::Point(map.ORIGIN_X, 0), cv::Point(map.ORIGIN_X, map.HEIGHT), cv::Scalar(80, 80, 80), 1, cv::LINE_8);
  for (auto p: XYs) {
    int ix = p.x / map.csize + map.ORIGIN_X;
    int iy =-p.y / map.csize + map.ORIGIN_Y;
    if (0 <= ix && ix < map.WIDTH && 0 <= iy && iy < map.HEIGHT) {
      cv::circle(img, cv::Point(ix, iy), 2, cv::Scalar(0, 0, 0), -1, cv::LINE_8);
    }
  }

  double min_dist = 99999;
  double best_x = 0.0;
  double best_y = 0.0;
  double best_a = 0.0;

  int separate = (xmax - xmin)/dd + 1;

  for (int ix = 0; ix < ((xmax - xmin)/dd + 1); ix++) {
    double rx = ix * dd + xmin;
    for (int iy = 0; iy < ((ymax - ymin)/dd + 1); iy++) {
      double ry = iy * dd + ymin;
      double min_dist_current_pos = 99999;
      double best_x_current_pos;
      double best_y_current_pos;
      double best_a_current_pos;
      for (double ra = 0.0; ra < 2*M_PI; ra += 5*2*M_PI/360.0) {
        double eval = evaluate(rx, ry, ra, data, skip, map, root, img);
        if (min_dist_current_pos > eval) {
          min_dist_current_pos = eval;
          best_a_current_pos = ra;
        }
      }
      match.emplace_back(rx, ry, best_a_current_pos, min_dist_current_pos);
      if (min_dist > min_dist_current_pos) {
        min_dist = min_dist_current_pos;
        best_x = rx;
        best_y = ry;
      }
    }
  }
  std::cout << best_x << " " << best_y << "\n";
  for (int i = 0; i < data.dat.size(); i++) {
    if (data.dat[i] > 35000) continue;
    double x = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
    double y = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

    // 座標変換
    double cs = cos(best_a);
    double sn = sin(best_a);
    double cx = cs * x - sn * y + best_x;
    double cy = sn * x + cs * y + best_y;
    int ax =  cx / map.csize + map.ORIGIN_X;
    int ay = -cy / map.csize + map.ORIGIN_Y;
    cv::circle(img, cv::Point(ax, ay), 2, cv::Scalar(0, 200, 0), -1, cv::LINE_8);
  }
  cv::imshow("BEST", img);
  //cv::waitKey();

  // 評価値の分布をファイルに出力
  int count = 0;
  std::ofstream fout("log");
  for (auto m: match) {
    count++;
    fout << m.x << " " << m.y << " " << m.match << "\n";
    if (count % separate == 0) fout << "\n";
  }
  return EXIT_SUCCESS;
}
