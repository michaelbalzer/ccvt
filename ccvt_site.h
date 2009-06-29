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

#ifndef CCVT_SITE_H
#define CCVT_SITE_H

#include <list>
#include <vector>

namespace ccvt {

  template<class Point>
  struct Site {

    typedef std::list<Site>     List;
    typedef std::vector<Site>   Vector;
    typedef std::vector<Site*>	VectorPtr;
    
    Site()
      : id(-1), 
        capacity(0) {
    }

    Site(const int id, const int capacity, const Point& location)
      : id(id), 
        capacity(capacity), 
        location(location) {
    }

    int	  id;
    int   capacity;
    Point location;

  };

}

#endif  // CCVT_SITE_H
