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

#ifndef CCVT_OPTIMIZER_H
#define CCVT_OPTIMIZER_H

#include <algorithm>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <math.h>

namespace ccvt {

  template<class Point, class Metric>
  class Optimizer {

  public:

    typedef Site<Point> Site;

  private:

    struct Entry;

    class KdTree {

    public:

      KdTree(std::vector<Entry> &entries, const Metric& metric, const int dimensions)
        : dimensions_(dimensions), root_(NULL), metric_(metric), lastNearestNeighborNode_(NULL) {
        if (!entries.empty()) {
          std::vector<Entry*> entriesPtr(entries.size());
          for (int i = 0; i < static_cast<int>(entriesPtr.size()); ++i) {
            entriesPtr[i] = &entries[i];
          }
          build(root_, entriesPtr, 0, static_cast<int>(entriesPtr.size() - 1), 0);
        }
      }
      
      ~KdTree() {
        if (root_ != NULL) {
          delete root_;
        }
      }

      bool empty() {
        return root_ == NULL || root_->disabled;
      }

      Entry* nearest_neighbor(const Point& point) {
        if (root_ == NULL || root_->disabled) {
          return NULL;
        }
        double minDistance = std::numeric_limits<double>::max();
        NearestNeighborResult result = nearest_neighbor(root_, point, 0, minDistance);      
        lastNearestNeighborNode_ = result.first;
        return result.second;
      }

      void disable_last_nearest_neighbor_node() {
        if (lastNearestNeighborNode_ != NULL) {
          lastNearestNeighborNode_->disabled = true;
          Node* parentNode = lastNearestNeighborNode_->parent;
          while (parentNode != NULL && parentNode->left->disabled && parentNode->right->disabled) {
              parentNode->disabled = true;
              parentNode = parentNode->parent;
          }
          lastNearestNeighborNode_ = NULL;
        }
      }

    private:

      KdTree(const KdTree& kdTree);
      KdTree& operator=(const KdTree& kdTree);

      struct Node {

        Node()
          : median(0), left(NULL), right(NULL), entry(NULL), parent(NULL), disabled(false) {
        }

        ~Node() {
          if (left != NULL) {
            delete left;
          }
          if (right != NULL) {
            delete right;
          }
        }

        double      median;
        Node*       left;
        Node*       right;
        Entry*      entry;
        Node*       parent;
        bool        disabled;

      };

      typedef std::pair<Node*, Entry*> NearestNeighborResult;

      struct Less {
        
        Less(int dimension) 
          : dimension(dimension) {
        }

        bool operator()(const Entry* p1, const Entry* p2) {
          return p1->site->location[dimension] < p2->site->location[dimension];
        }

        int dimension;

      };

      void build(Node*& node, std::vector<Entry*>& entries, const int min, const int max, const int depth) {
        node = new Node();      
        if (min == max) {
          node->entry = entries[min];
          return;
        }
        std::sort(entries.begin() + min, entries.begin() + max + 1, Less(depth % dimensions_));
        int medianIndex = (min + max) / 2;
        node->median = (entries[medianIndex]->site->location[depth % dimensions_] + entries[medianIndex + 1]->site->location[depth % dimensions_]) / 2;
        build(node->left, entries, min, medianIndex, depth + 1);
        build(node->right, entries, medianIndex + 1, max, depth + 1);
        node->left->parent = node;
        node->right->parent = node;
      }

      NearestNeighborResult nearest_neighbor(Node* node, const Point& point, const int depth, double& minDistance) {
        if (node->disabled) {
          return NearestNeighborResult(NULL, NULL);
        }

        if (node->entry != NULL) {
          minDistance = metric_.distance(point, node->entry->site->location);
          return NearestNeighborResult(node, node->entry);
        }

        Node* child;
        double otherDistance;
        if (point[depth % dimensions_] <= node->median) {
          child = node->left;
          otherDistance = node->median - point[depth % dimensions_];
        } else {
          child = node->right;
          otherDistance = point[depth % dimensions_] - node->median;
        }

        NearestNeighborResult result = nearest_neighbor(child, point, depth + 1, minDistance);

        if (minDistance > otherDistance) {
          double newMinDistance = minDistance;
          Node* otherChild = node->left;
          if (otherChild == child) {
            otherChild = node->right;
          }
          NearestNeighborResult newResult = nearest_neighbor(otherChild, point, depth + 1, newMinDistance);
          if (newMinDistance < minDistance) {
            result = newResult;
            minDistance = newMinDistance;
          }
        }

        return result;
      }

      const int     dimensions_;
      Node*         root_;
      const Metric& metric_;
      Node*         lastNearestNeighborNode_;

    };

  public:

    Optimizer(const Metric& metric)
      : metric_(metric), centroidal_(false) {
    }

    void clear() {
      entries_.clear();
      mapping_.clear();
    }

    void initialize(typename Site::Vector& sites, typename std::list<Point>& points, const bool centroidal) {
      centroidal_ = centroidal;
      int sitesSize = static_cast<int>(sites.size());
      entries_.reserve(sitesSize);

      for (int i = 0; i < sitesSize; ++i) {
        Site& site = sites[i];
        if (site.capacity > 0) {
          entries_.push_back(Entry(&site));
          Entry& entry = entries_.back();
          mapping_.insert(std::make_pair(site.id, &entry));
        } else {
          printf("\nWARNING: Site with zero capacity!");
        }
      }
      entries_.resize(entries_.size());

      KdTree kdTree(entries_, metric_, Point::D);
      while (!points.empty() && !kdTree.empty()) {
        const Point& point = points.back();
        Entry* entry = kdTree.nearest_neighbor(point);
        entry->points.push_back(point);
        if (static_cast<int>(entry->points.size()) == entry->site->capacity) {
          entry->points.resize(entry->points.size());
          kdTree.disable_last_nearest_neighbor_node();
        }
        points.pop_back();
      }

      for (unsigned int i = 0; i < entries_.size(); ++i) {
        entries_[i].energy = 0;
        int pointsSize = static_cast<int>(entries_[i].points.size());
        for (int j = 0; j < pointsSize; ++j) {
          entries_[i].energy += energy(entries_[i].points[j], entries_[i].site);
        }
        entries_[i].update(metric_);
      }

      if (!points.empty()) {
        printf("WARNING: Points not empty, %d remaining!", static_cast<int>(points.size()));
      }
    }

    bool optimize()
    {
      int entriesSize = static_cast<int>(entries_.size());
      std::vector<bool> stability(entriesSize, true);
      for (int i = 0; i < entriesSize; ++i) {
        for (int j = i + 1; j < entriesSize; ++j) {
          Entry* entry1 = &entries_[i];
          Entry* entry2 = &entries_[j];
          
          if (entry1->stable && entry2->stable ||
              metric_.distance(entry1->bounding.center, entry2->bounding.center) > entry1->bounding.radius + entry2->bounding.radius) {
            continue;
          }

          if (entry1->points.size() > entry2->points.size()) {
            std::swap(entry1, entry2);
          }

          Site* site1 = entry1->site;
          Site* site2 = entry2->site;
          std::vector<Point>* points1 = &entry1->points;
          std::vector<Point>* points2 = &entry2->points;

          double maxSquaredRadius = std::max(entry1->bounding.squaredRadius, entry2->bounding.squaredRadius);

          typename Candidate::Vector candidates1(points1->size());
          int size = static_cast<int>(points1->size());
          int count = 0;
          for (int k = 0; k < size; ++k) {
            Point& point = (*points1)[k];
            if (metric_.distance_square(point, entry2->bounding.center) <= maxSquaredRadius) {
              candidates1[count++] = Candidate(&point, energy(point, site1), energy(point, site2));
            }
          }
          if (count == 0) {
            continue;
          }
          candidates1.resize(count);
          std::make_heap(candidates1.begin(), candidates1.end());

          double minEnergy = -(candidates1.front().energySelf - candidates1.front().energyOther);
          typename Candidate::Vector candidates2(points2->size());
          size = static_cast<int>(points2->size());
          count = 0;
          for (int k = 0; k < size; ++k) {
            Point& point = (*points2)[k];
            if (metric_.distance_square(point, entry1->bounding.center) <= maxSquaredRadius) {
              double eSelf = energy(point, site2);
              double eOther = energy(point, site1);
              if (eSelf - eOther > minEnergy) {
                candidates2[count++] = Candidate(&point, eSelf, eOther);
              }
            }
          }
          if (count == 0) {
            continue;
          }
          candidates2.resize(count);
          std::make_heap(candidates2.begin(), candidates2.end());

          int maxSwaps = static_cast<int>(std::min(candidates1.size(), candidates2.size()));
          int swaps;
          for (swaps = 0; swaps < maxSwaps; ++swaps) {
            Candidate& candidate1 = candidates1.front();
            Candidate& candidate2 = candidates2.front();
            if (candidate1.energySelf - candidate1.energyOther + candidate2.energySelf - candidate2.energyOther <= 0) {
              break;
            }
            std::swap(*candidate1.point, *candidate2.point);
            entry1->energy += candidate2.energyOther - candidate1.energySelf;
            entry2->energy += candidate1.energyOther - candidate2.energySelf;
            std::pop_heap(candidates1.begin(), candidates1.end() - swaps);
            std::pop_heap(candidates2.begin(), candidates2.end() - swaps);
          }

          if (swaps > 0) {
            stability[i] = false;
            stability[j] = false;
            if (centroidal_) {
              entry1->site->location = metric_.centroid(entry1->site->location, entry1->points);
              entry2->site->location = metric_.centroid(entry2->site->location, entry2->points);
            }
            entry1->update(metric_);
            entry2->update(metric_);
          }
        }
      }

      bool stable = true;
      for (int i = 0; i < entriesSize; ++i) {
        entries_[i].stable = stability[i];
        stable &= stability[i];
      }
      return stable;
    }

    double energy() const {
      double e = 0;
      int entriesSize = static_cast<int>(entries_.size());
      for (int i = 0; i < entriesSize; ++i) {
        e += entries_[i].energy;
      }
      return e;
    }

    const std::vector<Point>* points(const int id) const {
      typename Entry::MapPtr::const_iterator it = mapping_.find(id);
      if (it == mapping_.end()) {
        return NULL;
      }
      return &it->second->points;
    }

    bool stable(const int id) const {
      typename Entry::MapPtr::const_iterator it = mapping_.find(id);
      return it == mapping_.end() || it->second->stable;
    }

  private:

    struct Bounding
    {

      Bounding()
        : radius(0), squaredRadius(0) {
      }

      Bounding(const Point& center, const double& radius)
        : center(center), radius(radius), squaredRadius(radius * radius) {
      }

      void update(const Point& site, const typename std::vector<Point>& points, const Metric& metric) {
        center = site;
        squaredRadius = 0;
        int pointsSize = static_cast<int>(points.size());
        for (int i = 0; i < pointsSize; ++i) {
          squaredRadius = std::max(squaredRadius, metric.distance_square(center, points[i]));
        }
        radius = sqrt(squaredRadius);
      }

      Point   center;
      double  radius;
      double  squaredRadius;

    };

    struct Entry 
    {

      typedef std::map<int, Entry*>	MapPtr;
      typedef	std::vector<Entry>    Vector;

      Entry()
        : site(NULL), stable(false) {
      }

      Entry(Site *const site)
        : bounding(site->location, 0), site(site), stable(false) {
      }

      void update(const Metric& metric) {
        stable = false;
        bounding.update(site->location, points, metric);
      }

      Bounding           bounding;
      std::vector<Point> points;
      Site*		           site;
      bool		           stable;
      double             energy;

    };

    struct Candidate
    {

      typedef std::vector<Candidate> Vector;

      Candidate()
        : point(NULL), energySelf(0), energyOther(0) {
      }

      Candidate(Point *const point, const double& energySelf, const double& energyOther)
        : point(point), energySelf(energySelf), energyOther(energyOther) {
      }

      inline bool operator<(const Candidate& candidate) const {
        return energySelf - energyOther < candidate.energySelf - candidate.energyOther;
      }

      Point*  point;
      double  energySelf;
      double  energyOther;

    };

    inline double energy(const Point& point, const Site *const site) const {
      return metric_.distance_square(point, site->location);
    }

    Metric                  metric_;
    typename Entry::Vector	entries_;
    typename Entry::MapPtr  mapping_;
    bool                    centroidal_;

  };

}

#endif	// CCVT_OPTIMIZER_H
