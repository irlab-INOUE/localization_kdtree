#include <vector>

struct Point {
  double x;
  double y;
  Point() {};
  Point(double x, double y) {
    this->x = x;
    this->y = y;
  }
};

struct Node {
  Point pt;
  Node *parent;
  Node *left;
  Node *right;

  Node() {
    parent = nullptr;
    left   = nullptr;
    right  = nullptr;
  };
};

Node *makeKDTree(std::vector<Point> &data, int start_index, int end_index, int depth, Node *parent_node)
{
  if (start_index >= end_index) return nullptr;

  // data をdepth でソート depth: 0 -> x軸 / 1 -> y軸 でそれぞれソート
  std::sort(
    data.begin() + start_index,
    data.begin() + end_index,
    [depth](const Point a, const Point b){return ((depth % 2) == 0) ? a.x < b.x : a.y < b.y;}
  );

  int mid_index = (start_index + end_index)/2;
  Node *nd   = new Node;
  nd->pt     = data[mid_index];
  nd->parent = parent_node;
  nd->left   = makeKDTree(data, start_index,   mid_index, depth + 1, nd);
  nd->right  = makeKDTree(data, mid_index + 1, end_index, depth + 1, nd);

  return nd;
}

Node *nearest_neighbor_search_KDTree(const Point &a, Node *nd, int depth) {
  auto dist = [](Point p1, Point p2) { return std::hypot(p1.x - p2.x, p1.y - p2.y); };

  double comp_a, comp_nd;
  if ( (depth % 2) == 0) { // X軸で比較
    comp_a  = a.x;
    comp_nd = nd->pt.x;
  } else {                 // Y軸で比較
    comp_a  = a.y;
    comp_nd = nd->pt.y;
  }

  if ((nd->left != nullptr && comp_a < comp_nd) || (nd->left != nullptr && nd->right == nullptr)) {
    Node *d = nearest_neighbor_search_KDTree(a, nd->left, depth + 1);
    if (dist(nd->pt, a) < dist(d->pt, a))
      d = nd;
    if (nd->right != nullptr && (comp_nd - comp_a)*(comp_nd - comp_a) < dist(d->pt, a)) {
      Node *d2 = nearest_neighbor_search_KDTree(a, nd->right, depth + 1);
      if (dist(d2->pt, a) < dist(d->pt, a))
        d = d2;
    }
    return d;
  }
  else if (nd->right != nullptr) {
    Node *d = nearest_neighbor_search_KDTree(a, nd->right, depth + 1);
    if (dist(nd->pt, a) < dist(d->pt, a))
      d = nd;
    if (nd->left != nullptr && (comp_nd - comp_a)*(comp_nd - comp_a) < dist(d->pt, a)) {
      Node *d2 = nearest_neighbor_search_KDTree(a, nd->left, depth + 1);
      if (dist(d2->pt, a) < dist(d->pt, a))
        d = d2;
    }
    return d;
  }
  else
    return nd;
}
