#ifndef BLISS_HEAP_HH
#define BLISS_HEAP_HH

/*
  Copyright (c) 2006-2011 Tommi Junttila
  Released under the GNU General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License version 3
  as published by the Free Software Foundation.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace bliss {

/** \internal
 * \brief A capacity bounded heap data structure.
 */

class Heap {
  unsigned int N;
  unsigned int n;
  unsigned int *array;
  void upheap(unsigned int k);
  void downheap(unsigned int k);

public:
  /**
   * Create a new heap.
   * init() must be called after this.
   */
  Heap() {
    array = 0;
    n = 0;
    N = 0;
  }
  ~Heap();

  /**
   * Initialize the heap to have the capacity to hold \e size elements.
   */
  void init(const unsigned int size);

  /**
   * Is the heap empty?
   * Time complexity is O(1).
   */
  bool is_empty() const { return (n == 0); }

  /**
   * Remove all the elements in the heap.
   * Time complexity is O(1).
   */
  void clear() { n = 0; }

  /**
   * Insert the element \a e in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  void insert(const unsigned int e);

  /**
   * Remove and return the smallest element in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  unsigned int remove();

  /**
   * Get the number of elements in the heap.
   */
  unsigned int size() const { return n; }
};

} // namespace bliss

#endif
