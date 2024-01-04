#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>

struct Point { // also used as a vector
  double x = 0;
  double y = 0;

  Point(double x, double y) : x(x), y(y) {}

  bool operator==(const Point& second) const {
    return (x - second.x) < 0.000000001 && (y - second.y) < 0.000000001;
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

  Point rotateVector(double angle) {
    angle *= M_PI / 180;
    Point rotated(x * cos(angle) - y * sin(angle),
      x * sin(angle) + y * cos(angle));
    return rotated;
  }

  void rotatePoint(const Point& center, double angle);

  void reflect(const Point& center);

  void scale(const Point& center, double coefficient);
};

// operations with vectors
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
double operator*(const Point& first, const Point& second) { // scalar
  return first.x * second.x + first.y * second.y;
}
Point operator/(const Point& point, double factor) {
  return Point(point.x / factor, point.y / factor);
}

Point Point::normalize() const {
  Point copy = *this;
  if (length() == 0) {
    std::cerr << "normalizing 0-vector\n";
  }
  return copy / length();
}

double myDistance(const Point& first, const Point& second) {
  double delta_x = first.x - second.x;
  double delta_y = first.y - second.y;
  return std::sqrt(delta_x * delta_x + delta_y * delta_y);
}

double Point::length() const {
  return myDistance(Point(0, 0), *this);
}

bool isParallel(Point vector1, Point vector2) {
  vector1 = vector1.normalize();
  vector2 = vector2.normalize();
  return vector1 == vector2 || vector1 == -vector2;
}

Point Point::projectionOn(Point vector) const {
  return vector.normalize() * (*this * vector) / vector.length();
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
  Point point;
  Point vector; // length = 1

public:
  Line(const Point& first, const Point& second) : point(first),
    vector((second - first).normalize()) {}

  Line(double coefficient, double shift) : point(Point(0, shift)),
    vector(Point(1, coefficient).normalize()) {}

  Line(const Point& point, double coefficient) : point(point),
    vector(Point(1, coefficient).normalize()) {}

  double asFunction(double x) const { // вертикаль => +inf если правее или равна, -inf если левее
    if (vector.x == 0) {
      return x >= point.x ? std::numeric_limits<double>::infinity() :
        -std::numeric_limits<double>::infinity();
    }
    return (x - point.x) / vector.x * vector.y + point.y;
  }

  bool operator==(const Line& second) {
    return (vector == second.vector || vector == -second.vector) && (point == second.point || isParallel(point - second.point, vector));
  }
  bool operator!=(const Line& second) {
    return !(*this == second);
  }

  Point direction_vector() const {
    return vector;
  }

  Point normal_vector() {
    return Point(vector.y, -vector.x).normalize();
  }

  double distance(const Point& start) const {
    double free_coefficient = -(vector.y * point.x - vector.x * point.y);
    return std::abs(vector.y * start.x - vector.x * start.y +
      free_coefficient) / vector.length();
  }

  Point reflect(const Point& to_reflect) const {
    Point normal = vector.leftNormal();
    if (distance(to_reflect + normal) > distance(to_reflect - normal)) {
      normal = -normal;
    }

    return to_reflect + normal * (2 * distance(to_reflect));
  }
};

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

template <typename T>
bool isSimular(const std::vector<T>& first, size_t start,
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
  return isSimular(first, start, second, step, 1);
}

void fillEdgeLengths(std::vector<double>& edge_lengths,
  const std::vector<Point>& vertices) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    edge_lengths[i] = myDistance(vertices[i], vertices[(i + 1) % vertices.size()]);
  }
}

void fillAngles(std::vector<double>& angles,
  const std::vector<Point>& vertices) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    Point direction1 = vertices[(i + 1) % vertices.size()] - vertices[i];
    Point direction2 = vertices[((i + vertices.size()) - 1) % vertices.size()]
      - vertices[i];
    angles[i] = direction1.normalize() * direction2.normalize(); // достаточно ли косинусов? надеюсь...
  }
}

bool isIntersect(Point ray_start, Point point1, Point point2) {
  if (point2.y < point1.y) { 
    Point copy = point1;
    point1 = point2;
    point2 = copy;
  }
  if (ray_start.y <= point1.y || ray_start.y > point2.y) {
    return false;
  }
  Point intersection = point1 + (point2 - point1) * ((ray_start.y - point1.y) /
    (point2.y - point1.y));
  return intersection.x >= ray_start.x;
}

class Polygon : public Shape {
protected:
  std::vector<Point> vertices;
public:
  Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}

  template <typename... Args>
  Polygon(const Args&... args) : vertices({ args... }) {}

  size_t verticesCount() const {
    return vertices.size();
  }

  std::vector<Point> getVertices() const {
    return vertices;
  }

  bool isConvex() const {
    for (size_t i = 0; i < vertices.size(); ++i) {
      const Point& next = vertices[(i + 1) % vertices.size()];
      const Line edge(vertices[i], next);

      const Point& next2 = vertices[(i + 2) % vertices.size()];
      bool sign = next2.y >= edge.asFunction(next2.x);
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

  double perimeter() const final {
    double perimeter = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
      perimeter += myDistance(vertices[i], vertices[(i + 1) % vertices.size()]);
    }
    return perimeter;
  }

  double area() const final {
    double area = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
      const Point& now = vertices[i];
      const Point& next = vertices[(i + 1) % vertices.size()];
      area += 1.0 / 2.0 * (now.x * next.y - next.x * now.y);
    }
    return area;
  }

  bool isEqualTo(const Shape& another) const final;
  double similarityCoefficient(const Shape& another) const;

  bool isCongruentTo(const Shape& another) const final {
    return similarityCoefficient(another) == 1.0;
  }

  bool isSimilarTo(const Shape& another) const final {
    return similarityCoefficient(another) != 0.0;
  }

  bool containsPoint(const Point& point) const final { // Метод трассировки луча (в сторону оси х), учёт числа пересечений
    size_t intersections = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
      intersections += isIntersect(point, vertices[i], vertices[(i + 1) % vertices.size()]);
    }
    return intersections % 2 == 1;
  }

  void rotate(const Point& center, double angle) final {
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i].rotatePoint(center, angle);
    }
  }
  void reflect(const Point& center) final {
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i].reflect(center);
    }
  }
  void reflect(const Line& axis) final {
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i] = axis.reflect(vertices[i]);
    }
  }
  void scale(const Point& center, double coefficient) final {
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i].scale(center, coefficient);
    }
  }
};

class Ellipse : public Shape {
private:
  Polygon toRhomb() const {
    Point direction1 = (focus2 - focus1);
    if (direction1.length() == 0) {
      direction1 = Point(1, 0);
    }
    else {
      direction1 = direction1.normalize();
    }
    Point direction2 = direction1.leftNormal();
    double half_axis2 = normal_half_axis();
    return Polygon(center() + direction1 * half_axis,
      center() + direction2 * half_axis2,
      center() - direction1 * half_axis,
      center() - direction2 * half_axis2);
  }

protected:
  Point focus1;
  Point focus2;
  double half_axis;

public:
  Ellipse(Point focus1, Point focus2, double length_sum)
    : focus1(focus1), focus2(focus2), half_axis(length_sum / 2) {};

  std::pair<Point, Point> fosuses() const {
    return std::make_pair(focus1, focus2);
  }

  double focal_length() const {
    return myDistance(focus1, focus2);
  }

  double normal_half_axis() const {
    double focal = focal_length() / 2;
    return std::sqrt(half_axis * half_axis - focal * focal);
  }

  double eccentricity() const {
    return (focal_length() / 2) / half_axis;
  }

  Point center() const {
    return Point((focus1.x + focus2.x) / 2, (focus1.y + focus2.y) / 2);
  }

  std::pair<Line, Line> directrices() const {
    double focal = focal_length() / 2;
    Point point1 = center() + (focus1 - focus2) / focal_length() *
      ((half_axis * half_axis) / focal);
    Point point2 = center() - (focus1 - focus2) / focal_length() *
      ((half_axis * half_axis) / focal);
    Point normal = Line(focus1, focus2).normal_vector();
    return std::make_pair(Line(point1, point1 + normal),
      Line(point2, point2 + normal));
  }

  double perimeter() const final {
    return 4*half_axis* std::comp_ellint_2(eccentricity());
  }

  double area() const final {
    return M_PI * half_axis * normal_half_axis();
  }

  bool isEqualTo(const Shape& another) const final;

  bool isCongruentTo(const Shape& another) const final;

  bool isSimilarTo(const Shape& another) const final;

  bool containsPoint(const Point& point) const final {
    return myDistance(point, focus1) + myDistance(point, focus2) < 2.0 * half_axis;
  }

  void rotate(const Point& center, double angle) final {
    focus1.rotatePoint(center, angle);
    focus2.rotatePoint(center, angle);
  }
  void reflect(const Point& center) final {
    focus1.reflect(center);
    focus2.reflect(center);
  }
  void reflect(const Line& axis) final {
    focus1 = axis.reflect(focus1);
    focus2 = axis.reflect(focus2);
  }
  void scale(const Point& center, double coefficient) final {
    focus1.scale(center, coefficient);
    focus2.scale(center, coefficient);
    Polygon rhomb = toRhomb();
    rhomb.scale(center, coefficient);
    half_axis = myDistance(rhomb.getVertices()[0], rhomb.getVertices()[2]) / 2;
  }
};

struct Circle : public Ellipse {
  Circle(Point center, double radius) : Ellipse(center, center, radius * 2) {};
  double radius() const {
    return half_axis;
  }
  Point center() {
    return focus1;
  }
};

bool Polygon::isEqualTo(const Shape& another) const {
  if (typeid(another) == typeid(Ellipse) || typeid(another) == typeid(Circle)) {
    return false;
  }
  const Polygon& second = dynamic_cast<const Polygon&>(another);
  if (vertices.size() != second.vertices.size()) {
    return false;
  }

  for (int step : {-1, 1}) {
    for (size_t start = 0; start < vertices.size(); ++start) {
      if (isEqual(vertices, start, second.vertices, step)) {
        return true;
      }
    }
  }

  return false;
}

double Polygon::similarityCoefficient(const Shape& another) const { // 0 => not similar
  if (typeid(another) == typeid(Ellipse) || typeid(another) == typeid(Circle)) {
    return 0.0;
  }
  const Polygon& second = dynamic_cast<const Polygon&>(another);
  if (vertices.size() != second.vertices.size()) {
    return 0.0;
  }

  std::vector<double> edge_lengths1(vertices.size());
  fillEdgeLengths(edge_lengths1, vertices);
  std::vector<double> edge_lengths2(vertices.size());
  fillEdgeLengths(edge_lengths2, vertices);

  std::vector<double> angles1(vertices.size());
  fillAngles(angles1, vertices);
  std::vector<double> angles2(vertices.size());
  fillAngles(angles2, vertices);

  for (int step : {-1, 1}) {
    for (size_t start = 0; start < vertices.size(); ++start) {
      if (isEqual(angles1, start, angles2, step) &&
        isSimular(edge_lengths1, start, edge_lengths2, step,
          edge_lengths2[start] / edge_lengths1[start])) {
        return edge_lengths2[start] / edge_lengths1[start];
      }
    }
  }

  return 0.0;
}

bool Ellipse::isEqualTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) && typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return (toRhomb() == second.toRhomb());
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) && typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return toRhomb().isCongruentTo(second.toRhomb());
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  if (typeid(another) != typeid(Ellipse) && typeid(another) != typeid(Circle)) {
    return false;
  }
  const Ellipse& second = dynamic_cast<const Ellipse&>(another);
  return toRhomb().isSimilarTo(second.toRhomb());
}

struct Rectangle : public Polygon {
  Rectangle(Point point1, Point point3, double ratio) {
    Point hypotenuse = (point3 - point1);
    Point height_base = point1 + (1 / (1 + ratio * ratio)) * hypotenuse;
    Point normal = hypotenuse.leftNormal();
    Point point2 = height_base + normal * std::sqrt((height_base - point1).length() *
      (point3 - height_base).length());
    Point point4 = point2;
    point4.reflect((point1 + point3) / 2);
    vertices = { point1, point2, point3, point4 };
  }

  Point center() const {
    return (vertices[0] + vertices[2]) / 2;
  }

  std::pair<Line, Line> diagonals() const {
    return std::make_pair(Line(vertices[0], vertices[2]),
      Line(vertices[1], vertices[3]));
  }
};

struct Triangle : public Polygon {
  using Polygon::Polygon;

  Circle circumscribedCircle() const {
    double radius = 1;
    for (size_t i = 0; i < 3; ++i) {
      radius *= myDistance(vertices[i], vertices[(i + 1) % 3]);
    }
    radius /= 4 * area();
    Point mid = (vertices[0] + vertices[1]) / 2;
    Point normal = (vertices[2] - mid).projectionOn((vertices[1] - vertices[0]).leftNormal()).normalize();
    double leg = myDistance(mid, vertices[0]);

    return Circle(mid + normal * std::sqrt(radius * radius - leg * leg),
      radius);
  }

  Circle inscribedCircle() const {
    Point direction_side1 = (vertices[1] - vertices[0]).normalize();
    Point direction_side2 = (vertices[2] - vertices[0]).normalize();
    Point direction_bisectrix = (direction_side1 + direction_side2).normalize();
    double radius = area() / (perimeter() / 2);

    return Circle(vertices[0] + direction_bisectrix * radius /
      std::sqrt((1 - direction_side1 * direction_side2) / 2),
      radius);
  }

  Point centroid() const {
    Point centroid(0, 0);
    for (size_t i = 0; i < 3; ++i) {
      centroid = centroid + vertices[i];
    }
    centroid = centroid / 3;
    return centroid;
  }

  Point orthocenter() const {
    Point orthocenter = vertices[0];
    orthocenter.reflect(circumscribedCircle().center());
    orthocenter.reflect((vertices[1] + vertices[2]) / 2);
    return orthocenter;
  }

  Line EulerLine() const {
    return Line(circumscribedCircle().center(), orthocenter());
  }

  Circle ninePointsCircle() const {
    return Triangle((vertices[0] + vertices[1]) / 2,
      (vertices[1] + vertices[2]) / 2,
      (vertices[2] + vertices[0]) / 2).circumscribedCircle();
  }
};

struct Square : public Rectangle {
  Square(Point point1, Point point2) : Rectangle(point1, point2, 1) {}

  Circle circumscribedCircle() {
    return Circle(center(), myDistance(vertices[0], vertices[2]) / 2);
  }

  Circle inscribedCircle() {
    Circle inscribed_circle = circumscribedCircle();
    inscribed_circle.scale(inscribed_circle.center(), 1 / std::sqrt(2));
    return inscribed_circle;
  }
};
