/*
  Copyright (C) 2009 Michael Balzer (michael.balzer@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <time.h>
#include "ccvt_metric.h"
#include "ccvt_optimizer.h"
#include "ccvt_point.h"
#include "ccvt_site.h"

using namespace ccvt;

// discrete space with constant density;
// the points form a regular grid
void constant_regular_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
  double n = sqrt(static_cast<double>(numberOfPoints));
  for (int x = 0; x < n; ++x) {
    for (int y = 0; y < n; ++y) {
      double dx = x / n * torusSize;
      double dy = y / n * torusSize;
      points.push_back(Point2(dx, dy));
    }
  }
}

// discrete space with constant density;
// the points are randomly distributed
void constant_random_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
  for (int i = 0; i < numberOfPoints; ++i) {
    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
    points.push_back(Point2(x, y));
  }
}

// discrete space with the density function e^(-20x^2-20y^2)+0.2sin^2(PIx)sin^2(PIy);
// the points are generated via rejection sampling
void nonconstant_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
  const double E = 2.718281828459;
  const double PI = 3.141592653590;
  while (points.size() < static_cast<unsigned int>(numberOfPoints)) {
    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
    double p = pow(E, -20.0 * x * x - 20.0 * y * y) + 0.2 * sin(PI * x) * sin(PI * x) * sin(PI * y) * sin(PI * y);
    double r = static_cast<double>(rand() % RAND_MAX) / RAND_MAX;
    if (p >= r) {
      points.push_back(Point2((x + 1) / 2 * torusSize, (y + 1) / 2 * torusSize));
    }
  }
}

// export sites to an EPS image
bool save_eps(const char* filename, const Site<Point2>::Vector& sites, const double width, const double height, const double radius) {
  std::ofstream stream(filename, std::ios::out);
  if (stream.bad()) {
    return false;
  }

  stream << "%!PS-Adobe EPSF-3.0\n";
  stream << "%%HiResBoundingBox: " << 0.0 << " " << 0.0 << " " << width << " " << height << "\n";
  stream << "%%BoundingBox: " << 0 << " " << 0 << " " << static_cast<int>(width) << " " << static_cast<int>(height) << "\n";
  stream << "\n";
  stream << "%% Sites: " << sites.size() << "\n";
  stream << "\n";
  stream << "/radius { " << radius << " } def\n";
  stream << "\n";
  stream << "/p { radius 0 360 arc closepath fill } def\n";
  stream << "\n";
  stream << "0 0 0 setrgbcolor\n";
  stream << "\n";
  for (unsigned int i = 0; i < sites.size(); ++i) {
    stream << sites[i].location.x << " " << sites[i].location.y << " p\n";
  }
  stream << "\n";
  stream << "showpage\n";

  stream.close();
  return true;
}

int main(int, char*[]) {
  const int     NUMBER_SITES      = 256;
  const int     NUMBER_POINTS     = 1024 * NUMBER_SITES;
  const double  TORUS_SIZE        = 1000;
  const bool    CONSTANT_DENSITY  = true;
  const bool    CENTROIDAL        = true;
  const bool    RESULT_PRINT      = false;
  const bool    RESULT_FILE       = true;
  const char*   RESULT_FILENAME   = "result.eps";
  const double  RESULT_RADIUS     = 5;
  
  typedef Optimizer<Site<Point2>, Point2, MetricToroidalEuclidean2> Optimizer;

  // intializing the underlying discrete space
  Point2::List points;
  if (CONSTANT_DENSITY) {
    constant_regular_density(points, NUMBER_POINTS, TORUS_SIZE);
  } else {
    nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE);
  }

  // initializing the Voronoi sites with equal capacity
  unsigned int overallCapacity = static_cast<int>(points.size());
  Site<Point2>::List sites;
  for (int i = 0; i < NUMBER_SITES; ++i) {
    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
    int capacity = overallCapacity / (NUMBER_SITES - i);
    overallCapacity -= capacity;
    sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
  }

  clock_t start = clock();

  // initializing the CCVT
  clock_t startInitialization = clock();
  printf("initialization...");
  Optimizer optimizer;
  MetricToroidalEuclidean2 metric(Point2(TORUS_SIZE, TORUS_SIZE));
  optimizer.initialize(sites, points, metric);
  printf("done\n");
  clock_t endInitialization = clock();

  // optimization
  int iteration = 0;
  bool stable;
  do {
    printf("iteration %d...", ++iteration);
    stable = optimizer.optimize(CENTROIDAL);
    printf("done\n");
  } while (!stable);
  
  clock_t end = clock();

  const Site<Point2>::Vector& result = optimizer.sites();
  
  // writing the Voronoi sites to console
  if (RESULT_PRINT) {
    printf("\nresult:\n");
    for (unsigned int i = 0; i < result.size(); ++i) {
      printf("site %d: %f, %f\n", result[i].id, result[i].location.x, result[i].location.y);
    }
  }

  printf("\ninitialization time: %.3f sec\n", static_cast<double>(endInitialization - startInitialization) / CLOCKS_PER_SEC);
  printf("computation time: %.3f sec\n", static_cast<double>(end - start) / CLOCKS_PER_SEC);

  // writing the Voronoi sites to EPS file
  if (RESULT_FILE) {
    if (save_eps(RESULT_FILENAME, result, TORUS_SIZE, TORUS_SIZE, RESULT_RADIUS)) {
      printf("\nresult saved in '%s'\n", RESULT_FILENAME);
    } else {
      printf("\nresult could not be saved in '%s'\n", RESULT_FILENAME);
    }
  }
  
  printf("\n");

  return 0;
}
