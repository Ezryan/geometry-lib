#include <cmath>
#include <iostream>
//#include <numbers>
#include <vector>

namespace Constant {
const double pi = 3.1415926535897932384626;
const double epsilon = 0.00005;
}


struct Point {
  double x = 0;
  double y = 0;

  Point() = default;
  Point(double x, double y);

  bool operator==(const Point& other) const;
  bool operator!=(const Point& other) const;
};


Point::Point(double x, double y): x(x), y(y) {}

bool Point::operator==(const Point &other) const {
  return std::abs(x - other.x) < Constant::epsilon &&
      std::abs(y - other.y) < Constant::epsilon;
}

bool Point::operator!=(const Point &other) const { return !(*this == other); }


class Line {
  friend class Polygon;
  friend class Ellipse;
 private:
  double angular_coef_;
  double shift_;

 public:
  Line(const Point& first, const Point& second);
  Line(double angular_coef, double shift);
  Line(const Point& dot, double angular_coef);

  bool operator==(const Line& other) const;
  bool operator!=(const Line& other) const;
};


Line::Line(const Point& first, const Point& second) {
  if (std::abs(first.x - second.x) < Constant::epsilon) {
    angular_coef_ = std::numeric_limits<double>::infinity();
    shift_ = first.x;
    return;
  }
  angular_coef_ = (first.y - second.y) / (first.x - second.x);
  shift_ = (first.x * second.y - first.y * second.x) / (first.x - second.x);
}

Line::Line(double angular_coef, double shift): angular_coef_(angular_coef), shift_(shift) {}

Line::Line(const Point& dot, double angular_coef) {
  if (std::isinf(angular_coef)) {
    angular_coef_ = std::numeric_limits<double>::infinity();
    shift_ = dot.x;
    return;
  }
  angular_coef_ = angular_coef;
  shift_ = dot.y - angular_coef * dot.x;
}

bool Line::operator==(const Line& other) const {
  return Point(angular_coef_, shift_) == Point(other.angular_coef_, other.shift_);
}

bool Line::operator!=(const Line &other) const { return !(*this == other); }


class Shape {
 public:
  virtual double perimeter() const = 0;
  virtual double area() const = 0;

  virtual bool isEqualTo(const Shape& other) const = 0;
  virtual bool isCongruentTo(const Shape& another) const = 0;
  virtual bool isSimilarTo(const Shape& another) const = 0;
  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;

  virtual ~Shape() = default;
};


bool operator==(const Shape& first, const Shape& second) {
  return first.isEqualTo(second);
}


class Polygon: public Shape {
 private:
  std::vector<Point> vertex_;

  static double VectorProduct(const Point& first, const Point& second);
  static double Deg2Rad(double angle);
  static bool isEqVec(const std::vector<Point>& first, const std::vector<Point>& second);
  static void CircularShift(std::vector<Point>& vec);
  static void Reverse(std::vector<Point>& vec);
  static double GetAng(const Point& first, const Point& second, const Point& third);

 public:
  explicit Polygon(const std::vector<Point>& vertex);
  template <class... T>
  explicit Polygon(T... vertex);
  size_t verticesCount() const;
  const std::vector<Point>& getVertices() const;
  bool isConvex() const;

  double perimeter() const final;
  double area() const final;

  bool isEqualTo(const Shape& other) const final;
  bool isCongruentTo(const Shape& another) const final;
  bool isSimilarTo(const Shape& another) const final;
  bool containsPoint(const Point& point) const final;

  void rotate(const Point& center, double angle) final;
  void reflect(const Point& center) final;
  void reflect(const Line& axis) final;
  void scale(const Point& center, double coefficient) final;
};


double Polygon::Deg2Rad(double angle) { return angle * Constant::pi / 180; }

double Polygon::VectorProduct(const Point& first, const Point& second) {
  return first.x * second.y - first.y * second.x;
}

bool Polygon::isEqVec(const std::vector<Point>& first, const std::vector<Point>& second) {
  for (size_t i = 0; i < first.size(); ++i) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

void Polygon::CircularShift(std::vector<Point>& vec) {
  Point buf = vec[0];
  for (size_t i = 0; i < vec.size() - 1; ++i) {
    vec[i] = vec[i + 1];
  }
  vec[vec.size() - 1] = buf;
}

void Polygon::Reverse(std::vector<Point>& vec) {
  for (size_t i = 0; i < vec.size() / 2; ++i) {
    std::swap(vec[i], vec[vec.size() - i - 1]);
  }
}

double Polygon::GetAng(const Point& first, const Point& second, const Point& third) {
  Point vec1 = {first.x - second.x, first.y - second.y};
  Point vec2 = {third.x - second.x, third.y - second.y};
  double cos_ang = (vec1.x * vec2.x + vec1.y * vec2.y) / hypot(vec1.x, vec1.y) / hypot(vec2.x, vec2.y);

  return acos(cos_ang);
}

Polygon::Polygon(const std::vector<Point>& vertex): vertex_(vertex) {}

template <class... T>
Polygon::Polygon(T... vertex): Polygon(std::vector<Point>({vertex...})) {}

size_t Polygon::verticesCount() const { return vertex_.size(); }

const std::vector<Point>& Polygon::getVertices() const { return vertex_; }

bool Polygon::isConvex() const {
  double vec_prod;
  for (size_t i = 0; i < vertex_.size() - 1; ++i) {
    Point current(vertex_[(i + 1) % vertex_.size()].x - vertex_[i % vertex_.size()].x,
                  vertex_[(i + 1) % vertex_.size()].y - vertex_[i % vertex_.size()].y);
    Point next(vertex_[(i + 2) % vertex_.size()].x - vertex_[(i + 1) % vertex_.size()].x,
               vertex_[(i + 2) % vertex_.size()].y - vertex_[(i + 1) % vertex_.size()].y);
    if (i != 0 && VectorProduct(current, next) * vec_prod < 0) {
      return false;
    }
    vec_prod = VectorProduct(current, next);
  }
  return true;
}

double Polygon::perimeter() const {
  double ret = 0;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    Point current(vertex_[(i + 1) % vertex_.size()].x - vertex_[i % vertex_.size()].x,
                  vertex_[(i + 1) % vertex_.size()].y - vertex_[i % vertex_.size()].y);
    ret += hypot(current.x, current.y);
  }
  return ret;
}

double Polygon::area() const {
  double ret = 0;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    ret += VectorProduct(vertex_[i], vertex_[(i + 1) % vertex_.size()]) / 2;
  }
  return std::abs(ret);
}

bool Polygon::isEqualTo(const Shape& other) const {
  const Polygon* polygon = dynamic_cast<const Polygon*>(&other);
  if (!polygon) {
    return false;
  }
  if (polygon->vertex_.size() != vertex_.size()) {
    return false;
  }

  std::vector<Point> copy = polygon->vertex_;
  for (size_t i = 0; i < copy.size(); ++i) {
    if (isEqVec(copy, vertex_)) {
      return true;
    }
    CircularShift(copy);
  }

  Reverse(copy);
  for (size_t i = 0; i < copy.size(); ++i) {
    if (isEqVec(copy, vertex_)) {
      return true;
    }
    CircularShift(copy);
  }

  return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
  if (!isSimilarTo(another)) {
    return false;
  }

  const Polygon* polygon = dynamic_cast<const Polygon*>(&another);

  std::vector<Point> mas1;
  std::vector<Point> mas2;
  Point buf;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    buf =
        {vertex_[i].x - vertex_[(i + 1) % vertex_.size()].x,
         vertex_[i].y - vertex_[(i + 1) % vertex_.size()].y};
    mas1.push_back({hypot(buf.x, buf.y), 0});
  }

  for (size_t i = 0; i < vertex_.size(); ++i) {
    buf =
        {polygon->vertex_[i].x - polygon->vertex_[(i + 1) % vertex_.size()].x,
         polygon->vertex_[i].y - polygon->vertex_[(i + 1) % vertex_.size()].y};
    mas2.push_back({hypot(buf.x, buf.y), 0});
  }

  return Polygon(mas1) == Polygon(mas2);
}

bool Polygon::isSimilarTo(const Shape& another) const {
  const Polygon* polygon = dynamic_cast<const Polygon*>(&another);

  if (!polygon) {
    return false;
  }

  if (polygon->vertex_.size() != vertex_.size()) {
    return false;
  }

  std::vector<Point> mas1;
  std::vector<Point> mas2;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    mas1.push_back({GetAng(vertex_[(i - 1 + vertex_.size()) % vertex_.size()],
                           vertex_[i],
                           vertex_[(i + 1) % vertex_.size()]), 0});
  }

  for (size_t i = 0; i < vertex_.size(); ++i) {
    mas2.push_back({GetAng(polygon->vertex_[(i - 1 + vertex_.size()) % vertex_.size()],
                           polygon->vertex_[i],
                           polygon->vertex_[(i + 1) % vertex_.size()]), 0});
  }

  return Polygon(mas1) == Polygon(mas2);
}

bool Polygon::containsPoint(const Point& point) const {
  int crossings = 0;
  double t;
  for (size_t i = 0; i < vertex_.size(); i++) {
    Point subsequentVertex = vertex_[(i + 1) % vertex_.size()];
    if (((vertex_[i].y <= point.y) && (subsequentVertex.y > point.y)) ||
        ((vertex_[i].y > point.y) && (subsequentVertex.y <= point.y))) {
      t = (point.y - vertex_[i].y) / (subsequentVertex.y - vertex_[i].y);
      if (point.x < vertex_[i].x + t * (subsequentVertex.x - vertex_[i].x)) {
        ++crossings;
      }
    }
  }

  return (crossings & 1) != 0;
}

void Polygon::rotate(const Point& center, double angle) {
  double buffer_x;
  double buffer_y;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    buffer_x =
        (vertex_[i].x - center.x) * cos(Deg2Rad(angle)) - (vertex_[i].y - center.y) * sin(Deg2Rad(angle)) + center.x;
    buffer_y =
        (vertex_[i].x - center.x) * sin(Deg2Rad(angle)) + (vertex_[i].y - center.y) * cos(Deg2Rad(angle)) + center.y;

    vertex_[i].x = buffer_x;
    vertex_[i].y = buffer_y;
  }
}

void Polygon::reflect(const Point& center) {
  for (size_t i = 0; i < vertex_.size(); ++i) {
    vertex_[i].x = 2 * center.x - vertex_[i].x;
    vertex_[i].y = 2 * center.y - vertex_[i].y;
  }
}

void Polygon::reflect(const Line& axis) {
  double x, y;
  Point buffer;
  for (size_t i = 0; i < vertex_.size(); ++i) {
    x = (vertex_[i].x + axis.angular_coef_ * (vertex_[i].y - axis.shift_))
        / (1 + axis.angular_coef_ * axis.angular_coef_);
    y = axis.angular_coef_ * x + axis.shift_;
    buffer = {2 * x - 2 * vertex_[i].x, 2 * y - 2 * vertex_[i].y};
    vertex_[i].x += buffer.x;
    vertex_[i].y += buffer.y;
  }
}

void Polygon::scale(const Point& center, double coefficient) {
  for (size_t i = 0; i < vertex_.size(); ++i) {
    vertex_[i].x = (vertex_[i].x - center.x) * coefficient + center.x;
    vertex_[i].y = (vertex_[i].y - center.y) * coefficient + center.y;
  }
}


class Ellipse: public Shape {
 private:
  std::pair<Point, Point> focus_;
  double summary_distance_;
  double dist_;
  double less_half_;
  double big_half_;

  void update();

 public:

  Ellipse(Point first, Point second, double dist);

  std::pair<Point,Point> focuses() const;
  std::pair<Line, Line> directrices() const;
  double eccentricity() const;
  Point center() const;

  double perimeter() const override;
  double area() const override;

  bool isEqualTo(const Shape& other) const final;
  bool isCongruentTo(const Shape& another) const final;
  bool isSimilarTo(const Shape& another) const final;
  bool containsPoint(const Point& point) const final;

  void rotate(const Point& center, double angle) override;
  void reflect(const Point& center) override;
  void reflect(const Line& axis) override;
  void scale(const Point& center, double coefficient) override;
};


void Ellipse::update() {
  dist_ = hypot(focus_.first.x - focus_.second.x, focus_.first.y - focus_.second.y);
  less_half_ = sqrt(summary_distance_ * summary_distance_ - dist_ * dist_) / 2;
  big_half_ = summary_distance_ / 2;
}

Ellipse::Ellipse(Point first, Point second, double distance): focus_(first, second), summary_distance_(distance) {
  dist_ = hypot(focus_.first.x - focus_.second.x, focus_.first.y - focus_.second.y);
  less_half_ = sqrt(summary_distance_ * summary_distance_ - dist_ * dist_) / 2;
  big_half_ = summary_distance_ / 2;
}

std::pair<Point,Point> Ellipse::focuses() const { return focus_; }

std::pair<Line, Line> Ellipse::directrices() const {
  Point cent = center();
  Point vec(2 * (focus_.first.x - cent.x) / dist_, 2 * (focus_.first.y - cent.y) / dist_);
  Line buf(focus_.first, focus_.second);

  cent.x += vec.x * big_half_ / eccentricity();
  cent.y += vec.y * big_half_ / eccentricity();

  Line right_dir(cent, -1.0 / buf.angular_coef_);

  cent.x -= 2 * vec.x * big_half_ / eccentricity();
  cent.y -= 2 * vec.y * big_half_ / eccentricity();

  Line left_dir(cent, -1.0 / buf.angular_coef_);

  return {left_dir, right_dir};
}

double Ellipse::eccentricity() const { return sqrt(1 - pow(less_half_ / big_half_, 2)); }

Point Ellipse::center() const {
  return {focus_.first.x / 2 + focus_.second.x / 2, focus_.first.y / 2 + focus_.second.y / 2};
}

double Ellipse::perimeter() const {
  return Constant::pi * (3 * (less_half_ + big_half_) -
      sqrt((3 * less_half_ + big_half_) * (less_half_ + 3 * big_half_)));
}

double Ellipse::area() const {
  return Constant::pi * less_half_ * big_half_;
}

bool Ellipse::isEqualTo(const Shape& another) const {
  const Ellipse* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (!ellipse) {
    return false;
  }
  return Polygon(focus_.first, focus_.second) == Polygon(ellipse->focus_.first, ellipse->focus_.second) &&
      summary_distance_ == ellipse->summary_distance_;
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  const Ellipse* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (!ellipse) {
    return false;
  }
  bool buf =
      Polygon(focus_.first, focus_.second).isCongruentTo(Polygon(ellipse->focus_.first, ellipse->focus_.second));
  return buf && ellipse->summary_distance_ == summary_distance_;
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  const Ellipse* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (!ellipse) {
    return false;
  }
  bool buf =
      Polygon(focus_.first, focus_.second).isSimilarTo(Polygon(ellipse->focus_.first, ellipse->focus_.second));
  return summary_distance_ / ellipse->summary_distance_ == dist_ / ellipse->dist_ && buf;
}

bool Ellipse::containsPoint(const Point& point) const {
  double dist = hypot(focus_.first.x - point.x, focus_.first.y - point.y) +
      hypot(focus_.second.x - point.x, focus_.second.y - point.y);
  return dist <= summary_distance_;
}

void Ellipse::rotate(const Point& center, double angle) {
  Polygon buf(focus_.first, focus_.second);
  buf.rotate(center, angle);
  focus_.first = buf.getVertices()[0];
  focus_.second = buf.getVertices()[1];
}

void Ellipse::reflect(const Point& center) {
  Polygon buf(focus_.first, focus_.second);
  buf.reflect(center);
  focus_.first = buf.getVertices()[0];
  focus_.second = buf.getVertices()[1];
}

void Ellipse::reflect(const Line& axis) {
  Polygon buf(focus_.first, focus_.second);
  buf.reflect(axis);
  focus_.first = buf.getVertices()[0];
  focus_.second = buf.getVertices()[1];
}

void Ellipse::scale(const Point& center, double coefficient) {
  Polygon buf(focus_.first, focus_.second);
  buf.scale(center, coefficient);
  focus_.first = buf.getVertices()[0];
  focus_.second = buf.getVertices()[1];
  summary_distance_ *= coefficient;
  update();
}


class Circle: public Ellipse {
 private:
  Point center_;
  double radius_;

 public:
  Circle(Point center, double radius);

  double radius() const;

  double perimeter() const override;
  double area() const override;

  void rotate(const Point& center, double angle) override;
  void reflect(const Point& center) override;
  void reflect(const Line& axis) override;
  void scale(const Point& center, double coefficient) override;
};


Circle::Circle(Point center, double radius): Ellipse(center, {center.x + Constant::epsilon, center.y}, 2 * radius),
                                             center_(center), radius_(radius) {}

double Circle::radius() const { return radius_; }

double Circle::perimeter() const { return 2 * Constant::pi * radius_; }

double Circle::area() const { return Constant::pi * radius_ * radius_; }

void Circle::rotate(const Point& center, double angle) {
  Polygon buffer(center_);
  buffer.rotate(center, angle);
  center_ = buffer.getVertices()[0];
  Ellipse::rotate(center, angle);
}

void Circle::reflect(const Point& center) {
  Polygon buffer(center_);
  buffer.reflect(center);
  center_ = buffer.getVertices()[0];
  Ellipse::reflect(center);
}

void Circle::reflect(const Line& axis) {
  Polygon buffer(center_);
  buffer.reflect(axis);
  center_ = buffer.getVertices()[0];
  Ellipse::reflect(axis);
}

void Circle::scale(const Point& center, double coefficient) {
  Polygon buffer(center_);
  buffer.scale(center, coefficient);
  center_ = buffer.getVertices()[0];
  radius_ *= coefficient;
  Ellipse::scale(center, coefficient);
}


class Rectangle: public Polygon {
 private:
  Point fisrt_;
  Point second_;
  double ratio_;

  static double Rad2Deg(double ang);

 public:
  Rectangle(Point first, Point second, double ratio);

  Point center() const;
  std::pair<Line, Line> diagonals() const;
};


double Rectangle::Rad2Deg(double ang) { return ang / Constant::pi * 180; }

Rectangle::Rectangle(Point first, Point second, double ratio):
    Polygon({first, {0, 0}, second, {0, 0}}),
    fisrt_(first), second_(second), ratio_(std::min(1 / ratio, ratio)) {
  Polygon buf(first, second);
  buf.rotate(center(), Rad2Deg(atan(-(2 * ratio_) / (1 - ratio_ * ratio_))));

  auto& vertex = const_cast<std::vector<Point>&>(Polygon::getVertices());
  vertex[1] = buf.getVertices()[0];
  vertex[3] = buf.getVertices()[1];
}

Point Rectangle::center() const { return {(fisrt_.x + second_.x) / 2, (fisrt_.y + second_.y) / 2}; }

std::pair<Line, Line> Rectangle::diagonals() const {
  Line buf1(fisrt_, second_);
  Polygon buf(fisrt_, second_);
  buf.rotate(center(), -(2 * ratio_) / (1 - ratio_ * ratio_));
  Line buf2(buf.getVertices()[0], buf.getVertices()[1]);
  return {buf1, buf2};
}


class Square: public Rectangle {
  double side_;

 public:
  Square(Point first, Point second);

  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
};


Square::Square(Point first, Point second): Rectangle(first, second, 1),
                                           side_(hypot(first.x - second.x, first.y - second.y) / sqrt(2)) {}

Circle Square::circumscribedCircle() const { return Circle(center(), side_ / sqrt(2)); }

Circle Square::inscribedCircle() const { return Circle(center(), side_ / 2); }


class Triangle: public Polygon {
 private:
  static Point middle(Point first, Point second);

 public:
  using Polygon::Polygon;

  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
  Point centroid() const;
  Point orthocenter() const;
  Line EulerLine() const;
  Circle ninePointsCircle() const;
};


Point Triangle::middle(Point first, Point second) {
  return {first.x / 2 + second.x / 2, first.y / 2 + second.y / 2};
}

Circle Triangle::circumscribedCircle() const {
  Point first = getVertices()[0];
  Point second = getVertices()[1];
  Point third = getVertices()[2];

  Point cent;

  cent.x = -0.5 * (first.y * (second.x * second.x - third.x * third.x + second.y * second.y - third.y * third.y) +
      second.y * (third.x * third.x - first.x * first.x + third.y * third.y - first.y * first.y) +
      third.y * (first.x * first.x - second.x * second.x + first.y * first.y - second.y * second.y)) /
      (first.x * (second.y - third.y) + second.x * (third.y - first.y) + third.x * (first.y - second.y));

  cent.y = 0.5 * (first.x * (second.x * second.x - third.x * third.x + second.y * second.y - third.y * third.y) +
      second.x * (third.x * third.x - first.x * first.x + third.y * third.y - first.y * first.y) +
      third.x * (first.x * first.x - second.x * second.x + first.y * first.y - second.y * second.y)) /
      (first.x * (second.y - third.y) + second.x * (third.y - first.y) + third.x * (first.y - second.y));

  return Circle(cent,
                hypot(getVertices()[0].x - getVertices()[1].x, getVertices()[0].y - getVertices()[1].y) *
                    hypot(getVertices()[1].x - getVertices()[2].x, getVertices()[1].y - getVertices()[2].y) *
                    hypot(getVertices()[2].x - getVertices()[0].x, getVertices()[2].y - getVertices()[0].y) /
                    4 / area());
}

Circle Triangle::inscribedCircle() const {
  double first = hypot(getVertices()[1].x - getVertices()[2].x, getVertices()[1].y - getVertices()[2].y);
  double second = hypot(getVertices()[2].x - getVertices()[0].x, getVertices()[2].y - getVertices()[0].y);
  double third = hypot(getVertices()[0].x - getVertices()[1].x, getVertices()[0].y - getVertices()[1].y);

  Point cent;
  cent.x = (getVertices()[0].x * first + second * getVertices()[1].x + third * getVertices()[2].x) /
      (first + second + third);

  cent.y = (getVertices()[0].y * first + second * getVertices()[1].y + third * getVertices()[2].y) /
      (first + second + third);

  return Circle(cent, 2 * area() / perimeter());
}

Point Triangle::centroid() const { return  {(getVertices()[0].x + getVertices()[1].x + getVertices()[2].x) / 3,
                                            (getVertices()[0].y + getVertices()[1].y + getVertices()[2].y) / 3};}

Point Triangle::orthocenter() const {
  Point buf_cent = circumscribedCircle().center();
  Point vec1 = {getVertices()[0].x - buf_cent.x, getVertices()[0].y - buf_cent.y};
  Point vec2 = {getVertices()[1].x - buf_cent.x, getVertices()[1].y - buf_cent.y};
  Point vec3 = {getVertices()[2].x - buf_cent.x, getVertices()[2].y - buf_cent.y};
  Point vec = {vec1.x + vec2.x + vec3.x, vec1.y + vec2.y + vec3.y};
  return {buf_cent.x + vec.x, buf_cent.y + vec.y};
}

Line Triangle::EulerLine() const { return {circumscribedCircle().center(), orthocenter()}; }

Circle Triangle::ninePointsCircle() const {
  Circle ret = circumscribedCircle();
  ret.scale(orthocenter(), 0.5);
  return ret;
}
