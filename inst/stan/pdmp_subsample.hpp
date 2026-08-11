#ifndef PDMP_SUBSAMPLE_HPP
#define PDMP_SUBSAMPLE_HPP

#include <ostream>
#include <vector>

#ifdef _WIN32
#define PDMP_EXPORT __declspec(dllexport)
#else
#define PDMP_EXPORT __attribute__((visibility("default")))
#endif

namespace pdmp_subsample {
thread_local std::vector<int> indices;
}

extern "C" {
PDMP_EXPORT void pdmp_set_subsample_indices(const int* idx, int size) {
  if (size < 0) return;
  pdmp_subsample::indices.assign(idx, idx + size);
}

PDMP_EXPORT void pdmp_clear_subsample_indices() {
  pdmp_subsample::indices.clear();
}

PDMP_EXPORT int pdmp_get_subsample_size() {
  return static_cast<int>(pdmp_subsample::indices.size());
}

// C-facing getter uses zero-based n and returns the stored zero-based index.
PDMP_EXPORT int pdmp_get_subsample_index(int n) {
  if (n < 0 || n >= static_cast<int>(pdmp_subsample::indices.size())) return -1;
  return pdmp_subsample::indices[static_cast<std::size_t>(n)];
}
}

// Stan external functions use one-based n and require one-based data indices.
inline int pdmp_get_subsample_size(std::ostream*) {
  return pdmp_get_subsample_size();
}

inline int pdmp_get_subsample_index(int n, std::ostream*) {
  return pdmp_get_subsample_index(n - 1) + 1;
}

#endif
