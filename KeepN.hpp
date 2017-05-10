/*
KeepN v1.5
Copyright(c) 2017 Steven Shave

Distributed under the MIT license

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once
#include <cfloat>
#include <list>
#include <utility>

template <class T> class KeepNAscending {
private:
  float cutoff = FLT_MAX;

public:
  unsigned int n;
  std::list<std::pair<T, float>> best;
  explicit KeepNAscending(const unsigned int ssize) { n = ssize; };
  void insert(const T &input, const float score) {
    if (score > cutoff)
      return;
    for (auto k = best.begin(); k != best.end(); k++) {
      if (score < std::get<1>(*k)) {
        best.insert(k, std::make_pair(input, score));
        if (best.size() > n)
          best.pop_back();
        if (best.size() == n)
          cutoff = std::get<1>(best.back());
        return;
      }
    }
    best.push_back(std::make_pair(input, score));
    if (best.size() > n)
      best.pop_back();
    if (best.size() == n)
      cutoff = std::get<1>(best.back());
  };

  void insert(const T &&input, const float score) {
    if (score > cutoff)
      return;
    for (auto k = best.begin(); k != best.end(); k++) {
      if (score < std::get<1>(*k)) {
        best.insert(k, std::make_pair(std::move(input), score));
        if (best.size() > n)
          best.pop_back();
        if (best.size() == n)
          cutoff = std::get<1>(best.back());
        return;
      }
    }
    best.push_back(std::make_pair(std::move(input), score));
    if (best.size() > n)
      best.pop_back();
    if (best.size() == n)
      cutoff = std::get<1>(best.back());
  };
};

template <class T> class KeepNDescending {
private:
  unsigned int sortsize;

public:
  float cutoff = -FLT_MAX;
  std::list<std::pair<T, float>> best;
  explicit KeepNDescending(const unsigned int ssize) { sortsize = ssize; };
  void insert(const T &input, const float score) {
    if (score < cutoff)
      return;
    for (auto k = best.begin(); k != best.end(); k++) {
      if (score > std::get<1>(*k)) {
        best.insert(k, std::make_pair(input, score));
        if (best.size() > sortsize)
          best.pop_back();
        if (best.size() == sortsize)
          cutoff = std::get<1>(best.back());
        return;
      }
    }
    best.push_back(std::make_pair(input, score));
    if (best.size() > sortsize)
      best.pop_back();
    if (best.size() == sortsize)
      cutoff = std::get<1>(best.back());
  };

  void insert(const T &&input, const float score) {
    if (score < cutoff)
      return;
    for (auto k = best.begin(); k != best.end(); k++) {
      if (score > std::get<1>(*k)) {
        best.insert(k, std::make_pair(std::move(input), score));
        if (best.size() > sortsize)
          best.pop_back();
        if (best.size() == sortsize)
          cutoff = std::get<1>(best.back());
        return;
      }
    }
    best.push_back(std::make_pair(std::move(input), score));
    if (best.size() > sortsize)
      best.pop_back();
    if (best.size() == sortsize)
      cutoff = std::get<1>(best.back());
  };
};