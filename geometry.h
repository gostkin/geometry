//
// Created by gosktin on 25.10.16.
//

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <utility>

namespace Geometry {
    const double EPS = 10e-9;
    const double PI = M_PI;

    template <typename T>
    inline T tabs(T object) {
        return object >= 0 ? object : -object;
    }

    template <typename T>
    inline T tpow(T object, unsigned degree) {
        if (degree == 2)
            return object * object;

        T s = object;
        for (unsigned i = 0; i < degree; ++i)
            s *= object;

        return s;
    }

///////////////////////////////////////////////////////////////////////////////
//    POINT
///////////////////////////////////////////////////////////////////////////////

    struct Point {
        double x;
        double y;

        Point() : x(0), y(0) {}

        Point(double x_, double y_) : x(x_), y(y_) {}

        inline bool operator==(const Point &right) const {
            return (tabs(x - right.x) < EPS) && (tabs(y - right.y) < EPS);
        }

        inline bool operator!=(const Point &right) const {
            return !(*this == right);
        }
    };

///////////////////////////////////////////////////////////////////////////////
//    VECTOR
///////////////////////////////////////////////////////////////////////////////

    struct Vector {
        Point start;
        Point stop;

        Vector(const Point &one, const Point &two) : start(one), stop(two) {}

        Vector(const Vector &v) {
            start = v.start;
            stop = v.stop;
        }

        inline bool operator==(const Vector &right) const {
            return (start == right.start) && (stop == right.stop);
        }

        inline bool operator!=(const Vector &right) const {
            return !(*this == right);
        }

        inline Vector getCanonical() const {
            return Vector(Point(0, 0), Point(stop.x - start.x, stop.y - start.y));
        }

        double operator[](const Vector &second) {
            Point canonicalFirst = getCanonical().stop;
            Point canonicalSecond = second.getCanonical().stop;

            return canonicalFirst.x * canonicalSecond.y - canonicalFirst.y * canonicalSecond.x;
        }

        double operator()(const Vector &second) {
            Point canonicalFirst = getCanonical().stop;
            Point canonicalSecond = second.getCanonical().stop;

            return canonicalFirst.x * canonicalSecond.x + canonicalFirst.y * canonicalSecond.y;
        }

        double length() const {
            return sqrt(tpow(stop.x - start.x, 2) + tpow(stop.y - start.y, 2));
        }

        Vector operator-() const {
            Vector temp(*this);
            Point tt = temp.stop;
            temp.stop = temp.start;
            temp.start = tt;

            return temp;
        }

        Vector operator+(const Vector &right) {
            Vector temp(*this);
            Vector r = right.getCanonical();
            temp.stop.x += r.stop.x;
            temp.stop.y += r.stop.y;

            return temp;
        }

        Vector operator-(const Vector &right) {
            Vector temp(*this);
            return temp + -right;
        }

        Vector operator*(double k) {
            if (stop == start)
                return *this;

            Vector temp(*this);
            temp.stop.x = start.x + k * (stop.x - start.x);
            temp.stop.y = start.y + k * (stop.y - start.y);

            return temp;
        }

        bool intersect(const Vector &v) {
            return intersect_(start.x, stop.x, v.start.x, v.stop.x)
                   && intersect_(start.y, stop.y, v.start.y, v.stop.y)
                   && area(start, stop, v.start) * area(start, stop, v.stop) <= 0
                   && area(v.start, v.stop, start) * area(v.start, v.stop, stop) <= 0;
        }

        Vector &normalize();

    private:
        inline double area(const Point &a, const Point &b, const Point &c) {
            return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        }

        inline bool intersect_(double a, double b, double c, double d) {
            if (a > b)
                std::swap(a, b);
            if (c > d)
                std::swap(c, d);
            return std::max(a, c) <= std::min(b, d);
        }
    };

///////////////////////////////////////////////////////////////////////////////
//    LINE
///////////////////////////////////////////////////////////////////////////////


    class Line {
    public:
        Line(const Point &first, const Point &second) {
            assert(first != second);
            if (tabs(first.y - second.y) < EPS) {
                xCoefficient_ = 0;
                yCoefficient_ = 1;
                xShift_ = 0;
                yShift_ = first.y;
            } else if (tabs(first.x - second.x) < EPS) {
                xCoefficient_ = 1;
                yCoefficient_ = 0;
                xShift_ = first.x;
                yShift_ = 0;
            } else {
                xCoefficient_ = (second.y - first.y) / (second.x - first.x);
                yCoefficient_ = 1;
                xShift_ = 0;
                yShift_ = first.y - xCoefficient_ * first.x;
            }
        }

        Line(double angularCoefficient, double shift) {
            xCoefficient_ = angularCoefficient;
            yCoefficient_ = 1;
            xShift_ = 0;
            yShift_ = shift;
        }

        Line(const Point &p, double angularCoefficient) {
            xCoefficient_ = angularCoefficient;
            yCoefficient_ = 1;
            xShift_ = 0;
            yShift_ = p.y - angularCoefficient * p.x;
        }

        inline bool operator==(const Line &right) const {
            return ((tabs(xCoefficient_ - right.xCoefficient_) < EPS) &&
                    (tabs(yCoefficient_ - right.yCoefficient_) < EPS) &&
                    (tabs(xShift_ - right.xShift_) < EPS) && (tabs(yShift_ - right.yShift_) < EPS));
        }

        inline bool operator!=(const Line &right) const {
            return !(*this == right);
        }

        Point intersectsAt(const Line &right) const {
            assert(right != *this);
            Point temp;
            if (right.yCoefficient_ == 0) {
                temp.x = right.xShift_;
                temp.y = yAt(temp.x);
            } else if (yCoefficient_ == 0) {
                temp.x = xShift_;
                temp.y = right.yAt(temp.x);
            } else {
                temp.x = (yShift_ - right.yShift_) / (right.xCoefficient_ - xCoefficient_);
                temp.y = yAt(temp.x);
            }
            return temp;
        }

        double yAt(double x) const {
            assert(yCoefficient_ != 0);
            return xCoefficient_ * x + yShift_;
        }

        double getAngularCoefficient() const {
            return xCoefficient_;
        }

        double getYCoefficient() const {
            return yCoefficient_;
        }

        double getShift() const {
            return yShift_;
        }

        double getXShift() const {
            return xShift_;
        }

        inline Line normal(Point point) const {
            if (getAngularCoefficient() != 0)
                return Line(point, (-1) / getAngularCoefficient());
            return Line(point, Point(point.x, point.y + 1));
        }

    private:
        double xCoefficient_;
        double yCoefficient_;
        double xShift_;
        double yShift_;
    };

///////////////////////////////////////////////////////////////////////////////
//    SHAPE
///////////////////////////////////////////////////////////////////////////////

    class Shape {
    public:
        virtual ~Shape() {}

        virtual double area() const = 0;
        virtual double perimeter() const = 0;

        virtual bool operator==(const Shape &right) const = 0;

        bool operator!=(const Shape &right) const {
            return !(*this == right);
        }

        virtual bool isCongruentTo(const Shape &right) const = 0;
        virtual bool isSimilarTo(const Shape &right) const = 0;
        virtual bool containsPoint(const Point &point) const = 0;

        virtual void rotate(const Point &center, double angle) = 0;
        virtual void reflex(const Point &center) = 0;
        virtual void reflex(const Line &axis) = 0;
        virtual void scale(const Point &center, double coefficient) = 0;
    };

///////////////////////////////////////////////////////////////////////////////
//    POLYGON
///////////////////////////////////////////////////////////////////////////////

    class Polygon : public Shape {
    public:
        template <typename ...Point>
        Polygon(const Point &... points) {
            points_ = {points...};
        }

        Polygon(const std::vector <Point> &points) : points_(points) {}

        inline size_t verticesCount() const {
            return points_.size();
        }

        std::vector <Point> getVertices() const {
            return points_;
        }

        bool isConvex() const {
            bool wasPositive = false, wasNegative = false;

            for (std::vector <Point>::const_iterator it = points_.begin(); it != points_.end() - 2; ++it)
                if (Vector(*it, *(it + 1))[Vector(*(it + 1), *(it + 2))] > 0)
                    wasPositive = true;
                else
                    wasNegative = true;

            if (Vector(*(points_.end() - 2), *(points_.end() - 1))[Vector(*(points_.end() - 1), *(points_.begin()))] > 0)
                wasPositive = true;
            else
                wasNegative = true;

            if (Vector(*(points_.end() - 1), *(points_.begin()))[Vector(*(points_.begin()), *(points_.begin() + 1))] > 0)
                wasPositive = true;
            else
                wasNegative = true;

            return !wasNegative == wasPositive;
        }

        double area() const {
            if (verticesCount() < 3)
                return 0;

            double sum = 0;
            size_t sz = points_.size();
            for (size_t it = 0; it != sz; ++it)
                sum += Vector(Point(0, 0), points_[it])[Vector(Point(0, 0), points_[(it + 1) % sz])] / 2;

            return tabs(sum);
        }

        double perimeter() const {
            double perimeter = 0;
            for (std::vector <Point>::const_iterator it = points_.begin(); it != points_.end() - 1; ++it)
                perimeter += Vector(*(points_.end() - 1), *(points_.begin())).length();

            perimeter += Vector(*(points_.end() - 1), *(points_.begin())).length();
            return perimeter;
        }

        bool operator==(const Shape &right) const {
            try {
                Polygon r = dynamic_cast<const Polygon &>(right);
                if (verticesCount() != r.verticesCount())
                    return false;

                if (verticesCount() == r.verticesCount() && verticesCount() == 3) {
                    return tabs(similarity(Polygon(r.points_[0], r.points_[1], r.points_[2])) - 1) < EPS ||
                           tabs(similarity(Polygon(r.points_[0], r.points_[2], r.points_[1])) - 1) < EPS ||
                           tabs(similarity(Polygon(r.points_[1], r.points_[0], r.points_[2])) - 1) < EPS;
                }
                size_t sz = points_.size();
                bool pep = false;
                for (size_t i = 0; i < sz; ++i) {
                    bool ok = true;
                    size_t iter = i;
                    for (size_t j = 0; j < sz; ++j) {
                        if (points_[iter % sz] == r.points_[j]) {
                            ++iter;
                        } else {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pep = true;
                        break;
                    }
                }

                std::reverse(r.points_.begin(), r.points_.end());

                for (size_t i = 0; i < sz; ++i) {
                    bool ok = true;
                    size_t iter = i;
                    for (size_t j = 0; j < sz; ++j) {
                        if (points_[iter % sz] == r.points_[j]) {
                            ++iter;
                        } else {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pep = true;
                        break;
                    }
                }

                return pep;
            } catch (...) {
                return false;
            }
        }

        bool isCongruentTo(const Shape &right) const {
            try {
                const Polygon &r = dynamic_cast<const Polygon &>(right);

                if (verticesCount() != r.verticesCount())
                    return false;

                if (verticesCount() == r.verticesCount() && verticesCount() == 3) {
                    return tabs(tabs(similarity(Polygon(r.points_[0], r.points_[1], r.points_[2]))) - 1) < EPS ||
                           tabs(tabs(similarity(Polygon(r.points_[0], r.points_[2], r.points_[1]))) - 1) < EPS ||
                           tabs(tabs(similarity(Polygon(r.points_[1], r.points_[0], r.points_[2])) - 1)) < EPS;
                }

                double t = similarity(r);
                return tabs(tabs(t) - 1) < EPS;
            } catch (...) {
                return false;
            }
        }

        bool isSimilarTo(const Shape &right) const {
            try {
                const Polygon &r = dynamic_cast<const Polygon &>(right);
                if (verticesCount() != r.verticesCount())
                    return false;

                if (verticesCount() == r.verticesCount() && verticesCount() == 3) {
                    return similarity(Polygon(r.points_[0], r.points_[1], r.points_[2])) ||
                           similarity(Polygon(r.points_[0], r.points_[2], r.points_[1])) ||
                           similarity(Polygon(r.points_[1], r.points_[0], r.points_[2]));
                }
                double t = similarity(r);

                return tabs(t) > EPS;
            } catch (...) {
                return false;
            }
        }

        bool containsPoint(const Point &point) const {
            Point p = points_[0], p1;
            size_t sz = points_.size();
            for (size_t i = 1; i < sz; ++i) {
                if (points_[i].x > p.x)
                    p = points_[i];
            }
            p1 = Point(p.x + 1, p.y + PI);

            size_t count = 0;
            for (size_t i = 0; i < sz; ++i) {
                if (Vector(points_[i], points_[(i + 1) % sz]).intersect(Vector(p1, point)))
                    ++count;
            }

            return count % 2 == 1;
        }

        void rotate(const Point &center, double angle) {
            angle *= PI / 180;
            for (std::vector <Point>::iterator point = points_.begin(); point != points_.end(); ++point) {
                point->y = (point->y - point->x * tan(angle)) / (tan(angle) * sin(angle) + cos(angle));
                point->x = point->x / cos(angle) + point->y * tan(angle);
            }
        }

        void reflex(const Point &center) {
            for (std::vector <Point>::iterator point = points_.begin(); point != points_.end(); ++point)
                *point = (Vector(Point(0, 0), center) - Vector(center, *point)).stop;
        }

        void reflex(const Line &axis) {
            for (std::vector <Point>::iterator point = points_.begin(); point != points_.end(); ++point) {
                Line l(*point, (-1) / axis.getAngularCoefficient());

                *point = (Vector(*point, l.intersectsAt(axis)) * 2).stop;
            }
        }

        void scale(const Point &center, double coefficient) {
            for (std::vector <Point>::iterator point = points_.begin(); point != points_.end(); ++point) {
                Vector temp = Vector(Point(0, 0), Point(point->x - center.x, point->y - center.y)) * coefficient;
                temp.stop.x += center.x;
                temp.stop.y += center.y;
                *point = temp.stop;
            }
        }

    protected:
        std::vector <Point> points_;


    private:
        double similarity(Polygon right) const {
            if (verticesCount() != right.verticesCount())
                return 0;

            for (size_t i = 0; i < right.verticesCount(); ++i) {
                double result = check(right);
                if (tabs(result - 0) > EPS)
                    return result;

                std::reverse(right.points_.begin(), right.points_.end());
                result = check(right);
                if (tabs(result - 0) > EPS)
                    return result;

                std::reverse(right.points_.begin(), right.points_.end());
                std::rotate(right.points_.begin(), right.points_.begin() + 1, right.points_.end());
            }

            return 0;
        }

        double check(Polygon right) const {
            bool equal = false;
            bool m = false;
            double coefficient =
                    Vector(right.points_[0], right.points_[1]).length() / Vector(points_[0], points_[1]).length();
            size_t sz = verticesCount();
            for (size_t i = 0; i < sz; ++i) {
                if (tabs(Vector(points_[i], points_[(i + 1) % sz]).length() * coefficient -
                         Vector(right.points_[i], right.points_[(i + 1) % sz]).length()) < EPS) {
                } else {
                    return 0;
                }

                double onev =
                        Vector(points_[i], points_[(i + 1) % sz])[Vector(points_[(i + 1) % sz],
                                                                         points_[(i + 2) % sz])] *
                        tpow(coefficient, 2);
                double twov = Vector(right.points_[i], right.points_[(i + 1) % sz])[Vector(right.points_[(i + 1) % sz],
                                                                                           right.points_[(i + 2) %
                                                                                                         sz])];
                double ones =
                        Vector(points_[i], points_[(i + 1) % sz])(
                                Vector(points_[(i + 1) % sz], points_[(i + 2) % sz])) *
                        tpow(coefficient, 2);
                double twos = Vector(right.points_[i], right.points_[(i + 1) % sz])(
                        Vector(right.points_[(i + 1) % sz], right.points_[(i + 2) % sz]));
                if (tabs(ones - twos) > EPS)
                    return 0;
                if (tabs(tabs(onev) - tabs(twov)) > EPS)
                    return 0;

                if (tabs(ones - twos) < EPS)
                    equal = true;
                else
                    m = true;
            }
            if (equal && m)
                return 0;

            if (m)
                return -coefficient;
            return coefficient;
        }
    };

    Vector &Vector::normalize() {
        Line l(Point(0, 0), start);
        Polygon p(stop);
        p.rotate(start, l.getAngularCoefficient() * 180 / PI);
        stop = p.getVertices()[0];
        return *this;
    }

///////////////////////////////////////////////////////////////////////////////
//    ELLIPSE
///////////////////////////////////////////////////////////////////////////////


    class Ellipse : public Shape {
    public:
        Ellipse(const Point &focusFirst = Point(0, 0), const Point &focusSecond = Point(0, 0), double distance = 0) :
                focusFirst_(focusFirst),
                focusSecond_(focusSecond),
                distance_(distance) {
        }

        Ellipse(double a, double b, Point center) {
            distance_ = 2 * a;
            focusFirst_ = focusSecond_ = center;
            focusFirst_.x -= sqrt(a * a - b * b);
            focusSecond_.x += sqrt(a * a - b * b);
        }

        std::pair <Point, Point> focuses() {
            return std::make_pair(focusFirst_, focusSecond_);
        };

        inline std::pair <Line, Line> directrixes() {
            Line l(focusFirst_, focusSecond_);
            Point one(center().x + getA() / eccentricity() * cos(atan(l.getAngularCoefficient())),
                      center().y + getA() / eccentricity() * sin(atan(l.getAngularCoefficient())));
            Point two(center().x - getA() / eccentricity() * cos(atan(l.getAngularCoefficient())),
                      center().y - getA() / eccentricity() * sin(atan(l.getAngularCoefficient())));

            return std::make_pair(l.normal(one), l.normal(two));
        };

        inline double eccentricity() {
            return getC() / getA();
        }

        Point center() const {
            return (Vector(focusFirst_, focusSecond_) * 0.5).stop;
        }

        double area() const {
            return PI * getA() * getB();
        }

        double perimeter() const {
            return PI * (3 * (getA() + getB()) - sqrt((3 * getA() + getB()) * (getA() + 3 * getB())));
        }

        bool operator==(const Shape &right) const {
            try {
                const Ellipse &r = dynamic_cast<const Ellipse &>(right);

                return tabs(getA() - r.getA()) < EPS && tabs(getB() - r.getB()) < EPS && center() == r.center();
            } catch (...) {
                return false;
            }
        }

        bool isCongruentTo(const Shape &right) const {
            try {
                const Ellipse &r = dynamic_cast<const Ellipse &>(right);
                return tabs(area() - r.area()) < EPS || tabs(getA() / r.getA() - getB() / r.getB()) < EPS ||
                       tabs(getA() / r.getB() - getB() / r.getA()) < EPS;
            } catch (...) {
                return false;
            }
        }

        bool isSimilarTo(const Shape &right) const {
            try {
                const Ellipse &r = dynamic_cast<const Ellipse &>(right);

                return tabs(getA() / r.getA() - getB() / r.getB()) < EPS;
            } catch (...) {
                return false;
            }
        }

        bool containsPoint(const Point &point) const {
            return Vector(point, focusFirst_).length() + Vector(point, focusSecond_).length() < distance_ + EPS;
        }

        void rotate(const Point &center, double angle) {
            angle *= PI / 180;

            focusFirst_.y = (focusFirst_.y - focusFirst_.x * tan(angle)) / (tan(angle) * sin(angle) + cos(angle));
            focusFirst_.x = focusFirst_.x / cos(angle) + focusFirst_.y * tan(angle);
            focusSecond_.y = (focusSecond_.y - focusSecond_.x * tan(angle)) / (tan(angle) * sin(angle) + cos(angle));
            focusSecond_.x = focusSecond_.x / cos(angle) + focusSecond_.y * tan(angle);
        }

        void reflex(const Point &center) {
            focusFirst_ = (Vector(Point(0, 0), center) - Vector(center, focusFirst_)).stop;
            focusSecond_ = (Vector(Point(0, 0), center) - Vector(center, focusSecond_)).stop;
        }

        void reflex(const Line &axis) {
            Line l(focusFirst_, (-1) / axis.getAngularCoefficient());

            focusFirst_ = (Vector(focusFirst_, l.intersectsAt(axis)) * 2).stop;
            focusSecond_ = (Vector(focusSecond_, l.intersectsAt(axis)) * 2).stop;
        }

        void scale(const Point &center, double coefficient) {
            Vector temp = Vector(Point(0, 0), Point(focusFirst_.x - center.x, focusFirst_.y - center.y)) * coefficient;
            temp.stop.x += center.x;
            temp.stop.y += center.y;
            focusFirst_ = temp.stop;

            temp = Vector(Point(0, 0), Point(focusSecond_.x - center.x, focusSecond_.y - center.y)) * coefficient;
            temp.stop.x += center.x;
            temp.stop.y += center.y;
            focusSecond_ = temp.stop;

            distance_ *= coefficient;
        }

        inline double getA() const {
            return distance_ / 2;
        }

        inline double getB() const {
            return sqrt(tpow(getA(), 2) - tpow(Vector(focusFirst_, focusSecond_).length() / 2, 2));
        }

        inline double getC() const {
            return sqrt(tpow(getA(), 2) - tpow(getB(), 2));
        }

    protected:
        Point focusFirst_;
        Point focusSecond_;
        double distance_;
    };


    class Circle : public Ellipse {
    public:
        Circle(Point center = Point(0, 0), double radius = 0) : Ellipse(center, center, radius * 2), radius_(radius) {
        }

        inline double radius() const {
            return radius_;
        }

        bool operator==(const Shape &right) const {

            try {
                const Circle &r = dynamic_cast<const Circle &>(right);
                return tabs(r.radius() - radius()) < EPS && center() == r.center();
            } catch (...) {
                try {
                    const Ellipse &r = dynamic_cast<const Ellipse &>(right);
                    return tabs(r.getA() - radius()) < EPS && getA() == getB();
                } catch (...) {
                    return false;
                }
            }
        }

        bool isCongruentTo(const Shape &right) const {
            try {
                const Circle &r = dynamic_cast<const Circle &>(right);
                return area() == r.area();
            } catch (...) {
                try {
                    const Ellipse &r = dynamic_cast<const Ellipse &>(right);
                    return tabs(r.getA() - radius()) < EPS && getA() == getB();
                } catch (...) {
                    return false;
                }
            }
        }

        bool isSimilarTo(const Shape &right) const {
            try {
                const Circle &r = dynamic_cast<const Circle &>(right);
                r.radius();
                return true;
            } catch (...) {
                try {
                    const Ellipse &r = dynamic_cast<const Ellipse &>(right);
                    return tabs(r.getA() - radius()) < EPS && getA() == getB();
                } catch (...) {
                    return false;
                }
            }
        }

        std::pair <Point, Point> intersectsAt(Line l) {
            double d = tpow(l.getAngularCoefficient(), 2) -
                       (tpow(l.getAngularCoefficient(), 2) + 1) * (tpow(l.getShift(), 2) - radius_ * radius());
            double x1 = (-l.getAngularCoefficient() + sqrt(d)) / (tpow(l.getAngularCoefficient(), 2) + 1);
            double x2 = (-l.getAngularCoefficient() - sqrt(d)) / (tpow(l.getAngularCoefficient(), 2) + 1);
            return std::make_pair(Point(x1, l.yAt(x1)), Point(x2, l.yAt(x2)));
        }

        bool containsPoint(const Point &point) const {
            return tabs(tpow(point.x - center().x, 2) + tpow(point.y - center().y, 2)) < tpow(radius_, 2);
        }

    private:
        double radius_;
    };


    class Rectangle : public Polygon {
    public:
        Rectangle(Point tl, Point br, double coefficient) {
            if (coefficient < 1)
                coefficient = 1 / coefficient;

            double k = tpow(Vector(tl, br).length(), 2);
            double one = k / (1 + tpow(coefficient, 2));
            double pointX = (br.x - tl.x) / k * tpow(one, 2);
            double pointY = (br.y - tl.y) / k * tpow(one, 2);
            pointX += (tl.y - br.y) / k * coefficient * tpow(one, 2);
            pointY += (tl.x - br.x) / k * coefficient * tpow(one, 2);
            points_.push_back(tl);
            points_.push_back(Point(tl.x + pointX, tl.y + pointY));
            points_.push_back(br);
            points_.push_back(Point(br.x - pointX, br.y - pointY));
        }

        Point center() const {
            Line l1(points_[0], points_[2]);
            Line l2(points_[1], points_[3]);
            return l1.intersectsAt(l2);
        }

        inline std::pair <Line, Line> diagonals() const {
            return std::make_pair(Line(points_[0], points_[2]), Line(points_[1], points_[3]));
        }
    };


    class Square : public Rectangle {
    public:
        Square(Point tl, Point br) : Rectangle(tl, br, 1) {}

        ~Square() {}

        Circle circumscribedCircle() const {
            std::pair <Line, Line> d = diagonals();
            return Circle(std::get <0>(d).intersectsAt(std::get <1>(d)), Vector(points_[0], points_[2]).length() / 2);
        }

        Circle inscribedCircle() const {
            std::pair <Line, Line> d = diagonals();
            return Circle(std::get <0>(d).intersectsAt(std::get <1>(d)), Vector(points_[0], points_[1]).length() / 2);
        }
    };


    class Triangle : public Polygon {
    public:
        Triangle(Point x, Point y, Point z) : Polygon(x, y, z) {}

        Point orthocenter() const {
            Point p = circumscribedCircle().center();
            return (Vector(p, points_[0]) + Vector(p, points_[1]) + Vector(p, points_[02])).stop;
        }

        Point centroid() const {
            return Point((points_[0].x + points_[1].x + points_[2].x) / 3,
                         (points_[0].y + points_[1].y + points_[2].y) / 3);
        }

        Circle circumscribedCircle() const {
            double d = 2 * (points_[0].x * (points_[1].y - points_[2].y) + points_[1].x * (points_[2].y - points_[0].y) +
                         points_[2].x * (points_[0].y - points_[1].y));
            double x = ((tpow(points_[0].x, 2) + tpow(points_[0].y, 2)) * (points_[1].y - points_[2].y) +
                         (tpow(points_[1].x, 2) + tpow(points_[1].y, 2)) * (points_[2].y - points_[0].y) +
                         (tpow(points_[2].x, 2) + tpow(points_[2].y, 2)) * (points_[0].y - points_[1].y)) / d;
            double y = ((tpow(points_[0].x, 2) + tpow(points_[0].y, 2)) * (points_[2].x - points_[1].x) +
                         (tpow(points_[1].x, 2) + tpow(points_[1].y, 2)) * (points_[0].x - points_[2].x) +
                         (tpow(points_[2].x, 2) + tpow(points_[2].y, 2)) * (points_[1].x - points_[0].x)) / d;

            return Circle(Point(x, y), Vector(points_[0], Point(x, y)).length());
        }

        Circle inscribedCircle() const {
            double a = Vector(points_[1], points_[2]).length(), b = Vector(points_[0], points_[2]).length(), c = Vector(
                    points_[1], points_[0]).length();
            double p = (a + b + c);

            return Circle(Point((a * points_[0].x + b * points_[1].x + c * points_[2].x) / p,
                                (a * points_[0].y + b * points_[1].y + c * points_[2].y) / p),
                          tabs(Vector(points_[0], points_[1])[Vector(points_[0], points_[2])]) / p);
        }

        inline Line EulerLine() const {
            return Line(circumscribedCircle().center(), orthocenter());
        }

        Circle ninePointsCircle() const {
            return Circle((Vector(circumscribedCircle().center(), orthocenter()) * 0.5).stop,
                          circumscribedCircle().radius() / 2);
        }
    };
}

#endif //GEOMETRY_GEOMETRY_H