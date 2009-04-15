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

#ifndef CCVT_FUNCTIONS_H
#define CCVT_FUNCTIONS_H

#include <fstream>
#include "ccvt_point.h"
#include "ccvt_site.h"

namespace ccvt {

  struct Functions {

    static bool save_eps(const char* filename, const Site<Point2>::Vector& sites, const double& width, const double& height, const double& radius) {
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

  };

}

#endif  // CCVT_FUNCTIONS_H
