#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>

namespace detail {
  double eps = 1e-9;
}

struct Point {
  double x = 0;
  double y = 0;

  Point(double x, double y) : x(x), y(y) {}

  bool operator==(const Point& second) const {
    return (x - second.x) < detail::eps && (y - second.y) < detail::eps;
  }
  bool operator!=(const Point& second) const {
    return !(*this == second);
  }

  double length() const;

  Point normalize() const;

  Point leftNormal() const {
    return Point(-y, x).normalize();
  }

  Point projectionOn(Point vector) const;

  Point rotateVector(double angle);

  void rotatePoint(const Point& center, double angle);

  void reflect(const Point& center);

  void scale(const Point& center, double coefficient);
};

Point Point::rotateVector(double angle) {
  angle *= M_PI / 180;
  Point rotated(x * cos(angle) - y * sin(angle),
    x * sin(angle) + y * cos(angle));
  return rotated;
}

namespace detail {
  double distance(const Point& first, const Point& second) {
    double delta_x = first.x - second.x;
    double delta_y = first.y - second.y;
    return std::sqrt(delta_x * delta_x + delta_y * delta_y);
  }
}

Point operator+(const Point& first, const Point& second) {
  return Point(first.x + second.x, first.y + second.y);
}
Point operator-(const Point& vector) {
  return Point(-vector.x, -vector.y);
}
Point operator-(const Point& first, const Point& second) {
  return Point(first.x - second.x, first.y - second.y);
}
Point operator*(const Point& point, double factor) {
  return Point(point.x * factor, point.y * factor);
}
Point operator*(double factor, const Point& point) {
  return Point(point.x * factor, point.y * factor);
}
double dotProduct(const Point& first, const Point& second) {
  return first.x * second.x + first.y * second.y;
}
Point operator/(const Point& point, double factor) {
  return Point(point.x / factor, point.y / factor);
}

Point Point::normalize() const {
  Point copy = *this;
  return copy / length();
}

double Point::length() const {
  return detail::distance(Point(0, 0), *this);
}

bool isParallel(Point first, Point second) {
  first = first.normalize();
  second = second.normalize();
  return first == second || first == -second;
}

Point Point::projectionOn(Point vector) const {
  return vector.normalize() * dotProduct(*this, vector) / vector.length();
}

void Point::rotatePoint(const Point& center, double angle) {
  *this = center + (*this - center).rotateVector(angle);
}

void Point::scale(const Point& center, double coefficient) {
  *this = center + coefficient * (*this - center);
}

void Point::reflect(const Point& center) {
  scale(center, -1);
}

class Line {
public:
  Line(const Point& first, const Point& second) : point(first),
    vector((second - first).normalize()) {}

  Line(double coefficient, double shift) : point(Point(0, shift)),
    vector(Point(1, coefficient).normalize()) {}

  Line(const Point& point, double coefficient) : point(point),
    vector(Point(1, coefficient).normalize()) {}

  bool operator==(const Line& second) {
    return (vector == second.vector || vector == -second.vector) &&
           (point == second.point || isParallel(point - second.point, vector));
  }

  bool operator!=(const Line& second) { return !(*this == second); }

  Point direction_vector() const { return vector; }

  Point normal_vector() { return Point(vector.y, -vector.x).normalize(); }

  double distance(const Point& start) const;

  Point reflect(const Point& to_reflect) const;

  double asFunction(double x) const;

private:
  Point point;
  Point vector;
};

double Line::distance(const Point& start) const {
  double free_coefficient = -(vector.y * point.x - vector.x * point.y);
  return std::abs(vector.y * start.x - vector.x * start.y +
    free_coefficient) / vector.length();
}

Point Line::reflect(const Point& to_reflect) const {
  Point normal = vector.leftNormal();
  if (distance(to_reflect + normal) > distance(to_reflect - normal)) {
    normal = -normal;
  }

  return to_reflect + normal * (2 * distance(to_reflect));
}

double Line::asFunction(double x) const {
  if (vector.x == 0) {
    return x >= point.x ? std::numeric_limits<double>::infinity() :
      -std::numeric_limits<double>::infinity();
  }
  return (x - point.x) / vector.x * vector.y + point.y;
}

class Shape {
public:
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool isEqualTo(const Shape&) const = 0;
  virtual bool isCongruentTo(const Shape&) const = 0;
  virtual bool isSimilarTo(const Shape&) const = 0;
  virtual bool containsPoint(const Point&) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;

  virtual ~Shape() = default;
};

bool operator==(const Shape& first, const Shape& second) {
  return first.isEqualTo(second);
}

namespace detail {
  template <typename T>
  bool isSimilar(const std::vector<T>& first, size_t start,
                 const std::vector<T>& second, int step, double ratio) {
    for (size_t i = 0; i < first.size(); ++i) {
      if (first[i] * ratio !=
        second[((start + second.size()) + i * step) % second.size()]) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool isEqual(const std::vector<T>& first, size_t start,
               const std::vector<T>& second, int step) {
    return detail::isSimilar(first, start, second, step, 1);
  }
}

void fillEdgeLengths(std::vector<double>& edge_lengths,
  const std::vector<Point>& vertices) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    edge_lengths[i] = detail::distance(vertices[i],
                                       vertices[(i + 1) % vertices.size()]);
  }
}

void fillAngles(std::vector<double>& angles,
  const std::vector<Point>& vertices) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    Point first_side = vertices[(i + 1) % vertices.size()] - vertices[i];
    Point second_side = vertices[((i + vertices.size()) - 1) % vertices.size()]
                        - vertices[i];
    angles[i] = dotProduct(first_side.normalize(), second_side.normalize());
  }
}

bool isIntersect(Point ray_start, Point first, Point second) {
  if (second.y < first.y) {
    Point copy = first;
    first = second;
    second = copy;
  }
  if (ray_start.y <= first.y || ray_start.y > second.y) {
    return false;
  }
  Point intersection = first + (second - first) * ((ray_start.y - first.y) /
                       (second.y - first.y));
  return intersection.x >= ray_start.x;
}

class Polygon : public Shape {
protected:
  std::vector<Point> vertices;
public:
  Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}

  template <typename... Args>
  Polygon(const Args&... args) : vertices({ args... }) {}

  size_t verticesCount() const { return vertices.size(); }

  std::vector<Point> getVertices() const { return vertices; }

  bool isConvex() const;

  double perimeter() const final;

  double area() const final;

  bool isEqualTo(const Shape& another) const final;
  double similarityCoefficient(const Shape& another) const;

  bool isCongruentTo(const Shape& another) const final {
    return similarityCoefficient(another) == 1.0;
  }

  bool isSimilarTo(const Shape& another) const final {
    return similarityCoefficient(another) != 0.0;
  }

  bool containsPoint(const Point& point) const final;

  void rotate(const Point& center, double angle) final;

  void reflect(const Point& center) final;

  void reflect(const Line& axis) final;

  void scale(const Point& center, double coefficient) final;
};

bool Polygon::isConvex() const {
  for (size_t i = 0; i < vertices.size(); ++i) {
    const Point& next = vertices[(i + 1) % vertices.size()];
    const Line edge(vertices[i], next);

    const Point& twice_next = vertices[(i + 2) % vertices.size()];
    bool sign = twice_next.y >= edge.asFunction(twice_next.x);
    for (size_t j = 0; j < vertices.size(); ++j) {
      if (j == i || j == (i + 1) % vertices.size()) {
        continue;
      }
      bool current_sign = vertices[j].y >= edge.asFunction(vertices[j].x);
      if (current_sign != sign) {
        return false;
      }
    }
  }
  return true;
}

double Polygon::perimeter() const {
  double perimeter = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    perimeter += detail::distance(vertices[i],
                                  vertices[(i + 1) % vertices.size()]);
  }
  return perimeter;
}

double Polygon::area() const {
  double area = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    const Point& now = vertices[i];
    const Point& next = vertices[(i + 1) % vertices.size()];
    area += 1.0 / 2.0 * (now.x * next.y - next.x * now.y);
  }
  return area;
}

bool Polygon::containsPoint(const Point& point) const {
  size_t intersections = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    intersections += isIntersect(point, vertices[i],
      vertices[(i + 1) % vertices.size()]);
  }
  return intersections % 2 == 1;
}

void Polygon::rotate(const Point& center, double angle) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].rotatePoint(center, angle);
  }
}
void Polygon::reflect(const Point& center) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].reflect(center);
  }
}
void Polygon::reflect(const Line& axis) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i] = axis.reflect(vertices[i]);
  }
}
void Polygon::scale(const Point& center, double coefficient) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].scale(center, coefficient);
  }
}

class Ellipse : public Shape {
public:
  Ellipse(Point first_focus, Point second_focus, double length_sum)
    : first_focus(first_focus), second_focus(second_focus),
      half_axis(length_sum / 2) {};

  std::pair<Point, Point> fosuses() const {
    return std::make_pair(first_focus, second_focus);
  }

  double focal_length() const {
    return detail::distance(first_focus, second_focus);
  }

  double normal_half_axis() const;

  double eccentricity() const { return (focal_length() / 2) / half_axis; }

  Point center() const {
    return Point((first_focus.x + second_focus.x) / 2,
                 (first_focus.y + second_focus.y) / 2);
  }

  std::pair<Line, Line> directrices() const;

  double perimeter() const final {
    return 4*half_axis* std::comp_ellint_2(eccentricity());
  }

  double area() const final { return M_PI * half_axis * normal_half_axis(); }

  bool isEqualTo(const Shape& another) const final;

  bool isCongruentTo(const Shape& another) const final;

  bool isSimilarTo(const Shape& another) const final;

  bool containsPoint(const Point& point) const final {
    return (detail::distance(point, first_focus) +
            detail::distance(point, second_focus)) < 2.0 * half_axis;
  }

  void rotate(const Point& center, double angle) final;

  void reflect(const Point& center) final;

  void reflect(const Line& axis) final;

  void scale(const Point& center, double coefficient) final;

protected:
  Point first_focus;
  Point second_focus;
  double half_axis;

private:
  Polygon toRhomb() const;
};

double Ellipse::normal_half_axis() const {
  double focal = focal_length() / 2;
  return std::sqrt(half_axis * half_axis - focal * focal);
}

std::pair<Line, Line> Ellipse::directrices() const {
  double focal = focal_length() / 2;
  Point point1 = center() + (first_focus - second_focus) / focal_length() *
    ((half_axis * half_axis) / focal);
  Point point2 = center() - (first_focus - second_focus) / focal_length() *
    ((half_axis * half_axis) / focal);
  Point normal = Line(first_focus, second_focus).normal_vector();
  return std::make_pair(Line(point1, point1 + normal),
    Line(point2, point2 + normal));
}

void Ellipse::rotate(const Point& center, double angle) {
  first_focus.rotatePoint(center, angle);
  second_focus.rotatePoint(center, angle);
}

void Ellipse::reflect(const Point& center) {
  first_focus.reflect(center);
  second_focus.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
  first_focus = axis.reflect(first_focus);
  second_focus = axis.reflect(second_focus);
}

void Ellipse::scale(const Point& center, double coefficient) {
  first_focus.scale(center, coefficient);
  second_focus.scale(center, coefficient);
  Polygon rhomb = toRhomb();
  rhomb.scale(center, coefficient);
  half_axis = detail::distance(rhomb.getVertices()[0],
    rhomb.getVertices()[2]) / 2;
}

Polygon Ellipse::toRhomb() const {
  Point direction = (second_focus - first_focus);
  if (direction.length() == 0) {
    direction = Point(1, 0);
  }
  else {
    direction = direction.normalize();
  }
  Point normal = direction.leftNormal();
  double second_half_axis = normal_half_axis();
  return Polygon(center() + direction * half_axis,
    center() + normal * second_half_axis,
    center() - direction * half_axis,
    center() - normal * second_half_axis);
}

struct Circle : public Ellipse {
  Circle(Point center, double radius) : Ellipse(center, center, radius * 2) {};
  double radius() const {
    return half_axis;
  }
  Point center() {
    return first_focus;
  }
};

bool Polygon::isEqualTo(const Shape& another) const {
  if (typeid(another) == typeid(Ellipse) ||
      typeid(another) == typeid(Circle)) {
    return false;
  }
  const Polygon& second = dynamic_cast<const Polygon&>(another);
  if (vertices.size() != second.vertices.size()) {
    return false;
  }

  for (int step : {-1, 1}) {
    for (size_t start = 0; start < vertices.size(); ++start) {
      if (detail::isEqual(vertices, start, second.vertices, step)) {
        return true;
      }
    }
  }

  return false;
}

double Polygon::similarityCoefficient(const Shape& another) const {
  if (typeid(another) == typeid(Ellipse) ||
      typeid(another) == typeid(Circle)) {
    return 0.0;
  }
  const Polygon& second = dynamic_cast<const Polygon&>(another);
  if (vertices.size() != second.vertices.size()) {
    return 0.0;
  }

  std::vector<double> edge_lengths_first(vertices.size());
  fillEdgeLengths(edge_lengths_first, vertices);
  std::vector<double> edge_lengths_second(vertices.size());
  fillEdgeLengths(edge_lengths_second, vertices);

  std::vector<double> angles_first(vertices.size());
  fillAngles(angles_first, vertices);
  std::vector<double> angles_second(vertices.size());
  fillAngles(angles_second, vertices);

  for (int step : {-1, 1}) {
    for (size_t start = 0; start < vertices.size(); ++start) {
      double ratio = edge_lengths_second[start] / edge_lengths_first[start];
      if (detail::isEqual(angles_first, start, angles_second, step) &&
          detail::isSimilar(edge_lengths_first, start,
                            edge_lengths_second, step, ratio)) {
        return edge_lengths_second[start] / edge_lengths_first[start];
      }
    }
  }

  return 0.0;
}

bool Ellipse::isEqualTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) &&
      typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return (toRhomb() == second.toRhomb());
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) &&
      typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return toRhomb().isCongruentTo(second.toRhomb());
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) &&
      typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return toRhomb().isSimilarTo(second.toRhomb());
}

struct Rectangle : public Polygon {
  Rectangle(Point first, Point third, double ratio);

  Point center() const { return (vertices[0] + vertices[2]) / 2; }

  std::pair<Line, Line> diagonals() const {
    return std::make_pair(Line(vertices[0], vertices[2]),
                          Line(vertices[1], vertices[3]));
  }
};

Rectangle::Rectangle(Point first, Point third, double ratio) {
  Point hypotenuse = (third - first);
  Point height_base = first + (1 / (1 + ratio * ratio)) * hypotenuse;
  Point normal = hypotenuse.leftNormal();
  Point second = height_base +
    normal * std::sqrt((height_base - first).length() *
      (third - height_base).length());
  Point fourth = second;
  fourth.reflect((first + third) / 2);
  vertices = { first, second, third, fourth };
}

struct Triangle : public Polygon {
  using Polygon::Polygon;

  Circle circumscribedCircle() const;

  Circle inscribedCircle() const;

  Point centroid() const;

  Point orthocenter() const;

  Line EulerLine() const {
    return Line(circumscribedCircle().center(), orthocenter());
  }

  Circle ninePointsCircle() const {
    return Triangle((vertices[0] + vertices[1]) / 2,
                    (vertices[1] + vertices[2]) / 2,
                    (vertices[2] + vertices[0]) / 2).circumscribedCircle();
  }
};

Circle Triangle::circumscribedCircle() const {
  double radius = 1;
  for (size_t i = 0; i < 3; ++i) {
    radius *= detail::distance(vertices[i], vertices[(i + 1) % 3]);
  }
  radius /= 4 * area();
  Point mid = (vertices[0] + vertices[1]) / 2;
  Point median = vertices[2] - mid;
  Point side = vertices[1] - vertices[0];
  Point normal = median.projectionOn(side.leftNormal()).normalize();
  double leg = detail::distance(mid, vertices[0]);

  return Circle(mid + normal * std::sqrt(radius * radius - leg * leg),
    radius);
}

Circle Triangle::inscribedCircle() const {
  Point first_direction = (vertices[1] - vertices[0]).normalize();
  Point secind_direction = (vertices[2] - vertices[0]).normalize();
  Point direction_bisectrix = (first_direction +
    secind_direction).normalize();
  double radius = area() / (perimeter() / 2);

  return Circle(vertices[0] + direction_bisectrix * radius /
    std::sqrt((1 - dotProduct(first_direction, secind_direction)) / 2),
    radius);
}

Point Triangle::centroid() const {
  Point centroid(0, 0);
  for (size_t i = 0; i < 3; ++i) {
    centroid = centroid + vertices[i];
  }
  centroid = centroid / 3;
  return centroid;
}

Point Triangle::orthocenter() const {
  Point orthocenter = vertices[0];
  orthocenter.reflect(circumscribedCircle().center());
  orthocenter.reflect((vertices[1] + vertices[2]) / 2);
  return orthocenter;
}

struct Square : public Rectangle {
  Square(Point first, Point second) : Rectangle(first, second, 1) {}

  Circle circumscribedCircle() {
    return Circle(center(), detail::distance(vertices[0], vertices[2]) / 2);
  }

  Circle inscribedCircle();
};

Circle Square::inscribedCircle() {
  Circle inscribed_circle = circumscribedCircle();
  inscribed_circle.scale(inscribed_circle.center(), 1 / std::sqrt(2));
  return inscribed_circle;
}
