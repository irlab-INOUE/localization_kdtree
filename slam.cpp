#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <opencv2/opencv.hpp>

#include "KDTree.h"
#include "Map.h"

#define DEBUG_DISPLAY
#define DELTA_D 1e-6

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

double evaluate(double rx, double ry, double ra, UrgData &data, int skip, Node *root) {
  double cs = cos(ra);
  double sn = sin(ra);
  double dist = 0;
  int count = 0;
  for (int i = 0; i < data.dat.size(); i+=skip) {
    if (data.dat[i] > 35000) continue;
    count++;
    double x = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
    double y = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

    // 座標変換
    double cx = cs * x - sn * y + rx;
    double cy = sn * x + cs * y + ry;
    Point addPoint(cx, cy);  // 探索する点

    // 最近傍点探索
    Node *nns = nearest_neighbor_search_KDTree(addPoint, root, 0);
    // 点間距離
    dist += std::hypot(cx - nns->pt.x, cy - nns->pt.y);
  }
  dist /= static_cast<double>(count);
  return dist;
}

struct GridMap {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double csize;

  int width;
  int height;

  int origin_x;
  int origin_y;

  std::vector<std::vector<std::vector<Point>>> cell;

  GridMap(double xmin, double xmax, double ymin, double ymax, double csize) {
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->csize = csize;

    width = (xmax - xmin)/csize;
    height = (ymax - ymin)/csize;
    cell.resize(height);
    for (int i = 0; i < height; i++) {
      cell[i].resize(width);
    }

    origin_x = abs(xmin)/csize;
    origin_y = abs(ymin)/csize;
  };

  void addPoint(double x, double y) {
    int ix = ( x - xmin)/csize;
    int iy = (-y - ymin)/csize;
    if (ix < 0 || ix > width || iy < 0 || iy > height) {
      ;
    } else {
      cell[iy][ix].emplace_back(x, y);
    }
  };

  void show() {
    cv::Mat img = cv::Mat(cv::Size(width, height), CV_8UC1, cv::Scalar(182, 182, 182));
    cv::line(img, cv::Point(0, origin_y), cv::Point(width, origin_y), cv::Scalar(200, 200, 200), 1, cv::LINE_8);
    cv::line(img, cv::Point(origin_x, 0), cv::Point(origin_x, height), cv::Scalar(200, 200, 200), 1, cv::LINE_8);

    for (int iy = 0; iy < height; iy++) {
      for (int ix = 0; ix < width; ix++) {
        if (cell[iy][ix].size() > 0) {
          img.at<uint8_t>(iy, ix) = 0;
        }
      }

    }
    cv::imshow("GridMap", img);
    cv::waitKey(5);
  };

  cv::Mat get_img() {
    cv::Mat img = cv::Mat(cv::Size(width, height), CV_8UC1, cv::Scalar(182, 182, 182));
    cv::line(img, cv::Point(0, origin_y), cv::Point(width, origin_y), cv::Scalar(200, 200, 200), 1, cv::LINE_8);
    cv::line(img, cv::Point(origin_x, 0), cv::Point(origin_x, height), cv::Scalar(200, 200, 200), 1, cv::LINE_8);

    for (int iy = 0; iy < height; iy++) {
      for (int ix = 0; ix < width; ix++) {
        if (cell[iy][ix].size() > 0) {
          img.at<uint8_t>(iy, ix) = 0;
        }
      }
    }
    return img;
  };
};

double fx(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (evaluate(x + dd, y, a, data, skip, root) - evaluate(x - dd, y, a, data, skip, root))/2.0/dd;
}
double fy(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (evaluate(x, y + dd, a, data, skip, root) - evaluate(x, y - dd, a, data, skip, root))/2.0/dd;
}
double fa(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (evaluate(x, y, a + dd, data, skip, root) - evaluate(x, y, a - dd, data, skip, root))/2.0/dd;
}
double fxx(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fx(x + dd, y, a, data, skip, root) - fx(x - dd, y, a, data, skip, root))/2.0/dd;
}
double fxy(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fx(x, y + dd, a, data, skip, root) - fx(x, y - dd, a, data, skip, root))/2.0/dd;
}
double fxa(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fx(x, y, a + dd, data, skip, root) - fx(x, y, a - dd, data, skip, root))/2.0/dd;
}
double fyy(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fy(x, y + dd, a, data, skip, root) - fy(x, y - dd, a, data, skip, root))/2.0/dd;
}
double fya(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fy(x, y, a + dd, data, skip, root) - fy(x, y, a - dd, data, skip, root))/2.0/dd;
}
double faa(double x, double y, double a, UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  return (fa(x, y, a + dd, data, skip, root) - fa(x, y, a - dd, data, skip, root))/2.0/dd;
}

std::tuple<double, double, double, double> descent_gradient(
  double current_x, double current_y, double current_a, 
  UrgData &data, int skip, Node *root) {
  double dd = DELTA_D;
  /* 勾配法 ---START--- */
  double eval = evaluate(current_x, current_y, current_a, data, skip, root);
  double tmpx = current_x;
  double tmpy = current_y;
  double tmpa = current_a;

  int loop = 0;
  while (1) {
    loop++;
    // 偏微分
    double dfx = fx(tmpx, tmpy, tmpa, data, skip, root);
    double dfy = fy(tmpx, tmpy, tmpa, data, skip, root);
    double dfa = fa(tmpx, tmpy, tmpa, data, skip, root);
    tmpx +=  - dfx * dd;
    tmpy +=  - dfy * dd;
    tmpa +=  - dfa * dd;
    if (std::hypot(dfx, dfy, dfa) < 1e-10 || loop > 10000) {
      eval = evaluate(tmpx, tmpy, tmpa, data, skip, root);
      current_x = tmpx;
      current_y = tmpy;
      current_a = tmpa;
      break;
    }
  }
  /* 勾配法---END--- */
  return std::tie(current_x, current_y, current_a, eval);
}

std::tuple<double, double, double, double> conjugate_gradient(
  double current_x, double current_y, double current_a, 
  UrgData &data, int skip, Node *root) {
	int itr_num = 100; 				// 繰り返し計算の上限
	double dfx, dfy, dfa; 			// 微分係数
	double p_dfx, p_dfy, p_dfa; 	// k-1回目の微分係数
	double mx, my, ma;
  double init_x = current_x;
  double init_y = current_y;
  double init_a = current_a;
	// 共役勾配法で極値を探索
	for (int k = 0; k < itr_num; k++) {
		// 偏微分係数を求める
		dfx = fx(current_x, current_y, current_a, data, skip, root);
		dfy = fy(current_x, current_y, current_a, data, skip, root);
		dfa = fa(current_x, current_y, current_a, data, skip, root);
		// 終了判定
		if ((dfx * dfx + dfy * dfy + dfa * dfa) < 1e-12) {
      double eval = evaluate(current_x, current_y, current_a, data, skip, root);
			return std::tie(current_x, current_y, current_a, eval);
		}
		// ヘッセ行列の成分と行列式
		double a11 = fxx(current_x, current_y, current_a, data, skip, root);
		double a12 = fxy(current_x, current_y, current_a, data, skip, root);
		double a13 = fxa(current_x, current_y, current_a, data, skip, root);
		double a21 = a12;
		double a22 = fyy(current_x, current_y, current_a, data, skip, root);
		double a23 = fya(current_x, current_y, current_a, data, skip, root);
		double a31 = a13;
		double a32 = a23;
		double a33 = faa(current_x, current_y, current_a, data, skip, root);
		double det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a12 * a21 * a33 - a11 * a23 * a32;
		if (std::isnan(det)) {
      current_x = init_x; current_y = init_y, current_a = init_a;
      double eval = evaluate(current_x, current_y, current_a, data, skip, root);
			return std::tie(current_x, current_y, current_a, eval);
		}
		double alpha;
		// 共役勾配法
		if (k == 0) {
			mx = dfx;
			my = dfy;
			ma = dfa;
		} else {
      // 定義通り
      alpha = -1.0 * (mx * (a11 * dfx + a12 * dfy + a13 * dfa) + my * (a21 * dfx + a22 * dfy + a23 * dfa)
        + ma * (a31 * dfx + a32 * dfy + a33 * dfa))/
        (mx * (a11 * mx + a12 * my + a13 * ma) + my * (a21 * mx + a22 * my + a23 * ma) + ma * (a31 * mx + a32 * my + a33 * ma));

			mx = dfx + alpha * mx;
			my = dfy + alpha * my;
			ma = dfa + alpha * ma;
		}
		p_dfx = dfx;
		p_dfy = dfy;
		p_dfa = dfa;

		// パラメータt
		// ここを直線探索にするとヘッセ行列は不要
		// これは2次近似
		double t = -1.0 * (mx * dfx + my * dfy + ma * dfa)
			/(mx * (a11 * mx + a12 * my + a13 * ma) + my * (a21 * mx + a22 * my + a23 * ma) + ma * (a31 * mx + a32 * my + a33 * ma));
		// 座標の更新
    current_x += t * mx;
    current_y += t * my;
    current_a += t * ma;
	}
  // 終了判定を失敗し繰り返し上限に達した場合
  current_x = init_x; current_y = init_y, current_a = init_a;
  double eval = evaluate(current_x, current_y, current_a, data, skip, root);
  return std::tie(current_x, current_y, current_a, eval);
}

void draw_measurement(cv::Mat &img, double csize, int origin_x, int origin_y, double x, double y, double a, UrgData &data) {
  double cs = cos(a);
  double sn = sin(a);
  for (int i = 0; i < data.dat.size(); i++) {
    if (data.dat[i] > 35000) continue;
    double lx = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
    double ly = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

    double cx = cs * lx - sn * ly + x;
    double cy = sn * lx + cs * ly + y;
    int ix = cx/csize + origin_x;
    int iy =-cy/csize + origin_y;
    cv::circle(img, cv::Point(ix, iy), 2, cv::Scalar(200, 200, 200), 1, cv::LINE_8);
  }
}

void draw_robot(cv::Mat &img, double csize, double origin_x, double origin_y, double x, double y, double a) {
  cv::circle(img, cv::Point(x/csize + origin_x, -y/csize + origin_y), 0.5/csize, cv::Scalar(20, 20, 20), 1, cv::LINE_AA);
  cv::line(img, 
           cv::Point(x/csize + origin_x, -y/csize + origin_y), 
           cv::Point(x/csize + origin_x + 0.5/csize*cos(a), -y/csize + origin_y - 0.5/csize*sin(a)),
           cv::Scalar(20, 20, 20), 1, cv::LINE_AA);
}

int main(int argc, char *argv[]) {
  GridMap gmap(-10.0, 10.0, -10.0, 10.0, 0.025);
  // 最初の計測を地図に登録
  std::ifstream fin("../urglog");

  double current_x = 0;
  double current_y = 0;
  double current_a = 0;
  double dd = 0.00001;
  const int skip = 10;        // LIDARデータのスキップ数

  std::string type;
  int count = 0;
  Node *root;   // KDTreeに変換した地図のルート
  double travel = 0.0;  // 走行距離, リセットする
  double total_travel = 0.0;  // 総走行距離 リセットしない
  double cumulative_rotate = 0.0;   // 地図の登録時点からの累積方向変化 リセットする
  // LIDARデータを逐次読み込み, ICP-SLAMを実行
  while (fin >> type) {
    UrgData data;
    data.type = type;
    fin 
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

    int start_index = 0;
    int last_index;
    Node *root_parent_is_nullptr = nullptr;
    Node *root;

    // 初回のLIDAR計測を地図登録
    if (count == 0) {
      // GridMapの作成
      for (int i = 0; i < data.dat.size(); i++) {
        if (data.dat[i] > 35000) continue;
        double x = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
        double y = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

        // 座標変換
        double cs = cos(current_a);
        double sn = sin(current_a);
        double cx = cs * x - sn * y + current_x;
        double cy = sn * x + cs * y + current_y;
        gmap.addPoint(cx, cy);
      }
      // 保持している地図からKDTreeを作成
      std::vector<Point> XYs;
      for (int y = 0; y < gmap.height; y++) {
        for (int x = 0; x < gmap.width; x++) {
          if (gmap.cell[y][x].size() > 0)
            XYs.emplace_back(gmap.cell[y][x][0].x, gmap.cell[y][x][0].y);   // 最初の登録点座標を用いるが，他にも考え方はありそう
            //XYs.emplace_back((x - gmap.origin_x) * gmap.csize, (-y + gmap.origin_y) * gmap.csize);
        }
      }
      last_index = XYs.size();
      root = makeKDTree(XYs, start_index, last_index, 0, root_parent_is_nullptr);
      count++;
      continue;
    }

    double pre_x = current_x;
    double pre_y = current_y;

    // SLAMの実行
    double eval;
#if 0
    // 最急降下法
    std::tie(current_x, current_y, current_a, eval)
      = descent_gradient(current_x, current_y, current_a, data, skip, root);
#endif
#if 1
    // 共役勾配法
    std::tie(current_x, current_y, current_a, eval)
      = conjugate_gradient(current_x, current_y, current_a, data, skip, root);
#endif
    double tt = std::hypot(current_x - pre_x, current_y - pre_y);
    travel += tt;
    total_travel += tt;
    std::cout << count << "\t" 
      << std::fixed
      << std::setprecision(3) << current_x << "\t" << current_y << "\t" 
      << std::setprecision(1) << current_a*180/M_PI <<  "\t" 
      << std::setprecision(3) << eval << " " 
      << (current_a - cumulative_rotate) * 180/M_PI << " "
      << travel << " " << total_travel << "\n";
    std::cout << std::defaultfloat;

    // 条件を満たせば地図へLIDAR計測の登録
    if (travel > 1.0 || abs(current_a - cumulative_rotate) > 10.0*M_PI/180) {
      travel = 0.0;
      cumulative_rotate = current_a;
      std::cout << "Updated Grid Map\n";
      // GridMapの作成
      for (int i = 0; i < data.dat.size(); i++) {
        if (data.dat[i] > 35000) continue;
        double x = data.dat[i] /1000.0 * cos((i * data.step + data.start_angle)*M_PI/180);
        double y = data.dat[i] /1000.0 * sin((i * data.step + data.start_angle)*M_PI/180);

        // 座標変換
        double cs = cos(current_a);
        double sn = sin(current_a);
        double cx = cs * x - sn * y + current_x;
        double cy = sn * x + cs * y + current_y;
        gmap.addPoint(cx, cy);
      }
      // 保持している地図からKDTreeを作成
      std::vector<Point> XYs;
      for (int y = 0; y < gmap.height; y++) {
        for (int x = 0; x < gmap.width; x++) {
          if (gmap.cell[y][x].size() > 0)
            XYs.emplace_back(gmap.cell[y][x][0].x, gmap.cell[y][x][0].y);   // 最初の登録点座標を用いるが，他にも考え方はありそう
            //XYs.emplace_back((x - gmap.origin_x) * gmap.csize, (-y + gmap.origin_y) * gmap.csize);
        }
      }
      last_index = XYs.size();
      root = makeKDTree(XYs, start_index, last_index, 0, root_parent_is_nullptr);
    }

    // 経過の表示
    cv::Mat img = gmap.get_img();
    draw_measurement(img, gmap.csize, gmap.origin_x, gmap.origin_y, current_x, current_y, current_a, data);  // 点群を座標変換して描画
    draw_robot(img, gmap.csize, gmap.origin_x, gmap.origin_y, current_x, current_y, current_a);              // ロボットの描画

    count++;
    //gmap.show();
    cv::imshow("Progress", img);
    cv::waitKey(5);
  }
  cv::imwrite("gmap.png", gmap.get_img());

  std::cout << "END" << std::endl;
  return EXIT_SUCCESS;
}
