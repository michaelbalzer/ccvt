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

#include <time.h>
#include "ccvt_functions.h"
#include "ccvt_metric.h"
#include "ccvt_optimizer.h"
#include "ccvt_point.h"
#include "ccvt_site.h"

using namespace ccvt;

int main(int , char * [])
{
  const int     NUMBER_SITES  = 1024;
  const int     NUMBER_POINTS = 1024;
  const double  TORUS_SIZE    = 1000;
  const bool    CENTROIDAL    = true;
  const bool    RESULT_PRINT  = false;
  const char*   RESULT_FILE   = "result.eps";
  const double  RESULT_RADIUS = 3;
  
  typedef Optimizer<Point2, MetricToroidalEuclidean2> Optimizer;

  Point2::List points;
  for (double x = 0; x < NUMBER_POINTS; ++x) {
    for (double y = 0; y < NUMBER_POINTS; ++y) {
      double dx = x / NUMBER_POINTS * TORUS_SIZE;
      double dy = y / NUMBER_POINTS * TORUS_SIZE;
      points.push_back(Point2(dx, dy));
    }
  }

  unsigned int overallCapacity = static_cast<int>(points.size());
  Optimizer::Site::Vector sites(NUMBER_SITES);
  for (int i = 0; i < static_cast<int>(sites.size()); ++i) {
    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
    int capacity = overallCapacity / (sites.size() - i);
    overallCapacity -= capacity;
    sites[i] = Optimizer::Site(i, capacity, Point2(x, y));
  }

  clock_t start = clock();

  clock_t startInitialization = clock();
  printf("initialization...");
  MetricToroidalEuclidean2 metric(Point2(TORUS_SIZE, TORUS_SIZE));
  Optimizer optimizer(metric);
  optimizer.initialize(sites, points, CENTROIDAL);
  printf("done\n");
  clock_t endInitialization = clock();

  int iteration = 0;
  bool stable;
  do {
    printf("iteration %d...", ++iteration);
    stable = optimizer.optimize();
    printf("done\n");
  } while (!stable);
  
  clock_t end = clock();
  
  if (RESULT_PRINT) {
    printf("\nresult:\n");
    for (unsigned int j = 0; j < sites.size(); ++j) {
      printf("site %d: %f, %f\n", sites[j].id, sites[j].location.x, sites[j].location.y);
    }
  }

  printf("\ninitialization time: %.3f sec\n", static_cast<double>(endInitialization - startInitialization) / CLOCKS_PER_SEC);
  printf("computation time: %.3f sec\n", static_cast<double>(end - start) / CLOCKS_PER_SEC);

  if (Functions::save_eps(RESULT_FILE, sites, TORUS_SIZE, TORUS_SIZE, RESULT_RADIUS)) {
    printf("\nresult saved in '%s'\n", RESULT_FILE);
  } else {
    printf("\nresult could not be saved in '%s'\n", RESULT_FILE);
  }
  
  printf("\n");

  return 0;
}
