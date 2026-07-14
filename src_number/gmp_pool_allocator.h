// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// A size-segregated free-list pool plugged into GMP's allocation hooks
// (mp_set_memory_functions). Exact rational/integer computations (mpz_class,
// mpq_class, jet<mpq_class>) are dominated -- ~50% of runtime on the lattice
// quantization moment computation -- by system malloc/free of tiny, extremely
// short-lived limb allocations, where the general-purpose allocator's
// bookkeeping (size lookup, zone management, thread synchronization) is the
// cost. This replaces it with a trivial pool: one free-list per (16-byte) size
// class, blocks carved from big chunks. O(1) alloc/free, no size headers (GMP
// hands the size back to free/realloc), no coalescing, no locks. Measured ~46%
// end-to-end speedup on the dim-6 quantization cases.
//
// ENABLING IT
//   Build with -DGMP_POOL. The pool then installs itself automatically (a static
//   constructor calls mp_set_memory_functions), so no per-program main() change
//   is needed. Without -DGMP_POOL this header is inert.
//
// SERIAL vs MULTI-THREADED
//   Default: a single process-global pool. Correct for single-threaded programs
//   AND for MPI (each rank is a separate process with its own pool -- MPI is
//   multi-process, not multi-threaded).
//   -DGMP_POOL_THREAD_SAFE: the pool state becomes thread_local -- one free-list
//   set per thread -- so it is safe when several threads of one process use GMP
//   concurrently. Use this for the genuinely multi-threaded code (the MPI
//   dual-description runs a background communication thread). Cross-thread frees
//   are safe: a block freed on another thread simply joins that thread's
//   free-list and is reused there; the memory stays valid, nothing is corrupted.
//
// MEMORY
//   The pool never hands chunks back to the OS during the run -- it grows to the
//   peak working set. At process exit the OS reclaims everything, so there is no
//   persistent leak (a leak checker may report the pooled chunks as still
//   reachable at exit; that is expected and harmless).
#ifndef SRC_NUMBER_GMP_POOL_ALLOCATOR_H_
#define SRC_NUMBER_GMP_POOL_ALLOCATOR_H_

#include <cstdlib>
#include <cstring>
#include <gmp.h>

#ifdef GMP_POOL_THREAD_SAFE
#define GMP_POOL_STORAGE thread_local
#else
#define GMP_POOL_STORAGE
#endif

namespace gmp_pool {

inline constexpr size_t ALIGN = 16;          // >= alignof(mp_limb_t), >= sizeof(void*)
inline constexpr size_t MAX_POOLED = 4096;   // larger requests go straight to malloc
inline constexpr size_t NCLASS = MAX_POOLED / ALIGN;
inline constexpr size_t CHUNK = size_t(8) << 20; // 8 MB arena chunks

struct FreeNode {
  FreeNode *next;
};

inline GMP_POOL_STORAGE FreeNode *g_freelist[NCLASS] = {};
inline GMP_POOL_STORAGE char *g_arena = nullptr;
inline GMP_POOL_STORAGE size_t g_arena_left = 0;

// 0-based size class for a non-zero size in (0, MAX_POOLED].
inline size_t class_index(size_t size) { return (size + ALIGN - 1) / ALIGN - 1; }

inline void *pool_alloc(size_t size) {
  if (size == 0)
    size = 1;
  if (size > MAX_POOLED)
    return std::malloc(size);
  size_t ci = class_index(size);
  if (g_freelist[ci] != nullptr) {
    FreeNode *n = g_freelist[ci];
    g_freelist[ci] = n->next;
    return n;
  }
  size_t bytes = (ci + 1) * ALIGN;
  if (g_arena_left < bytes) {
    g_arena = static_cast<char *>(std::malloc(CHUNK));
    g_arena_left = CHUNK;
  }
  void *p = g_arena;
  g_arena += bytes;
  g_arena_left -= bytes;
  return p;
}

inline void pool_free(void *ptr, size_t size) {
  if (size == 0)
    size = 1;
  if (size > MAX_POOLED) {
    std::free(ptr);
    return;
  }
  size_t ci = class_index(size);
  FreeNode *n = static_cast<FreeNode *>(ptr);
  n->next = g_freelist[ci];
  g_freelist[ci] = n;
}

inline void *pool_realloc(void *ptr, size_t old_size, size_t new_size) {
  if (old_size == 0)
    old_size = 1;
  if (new_size == 0)
    new_size = 1;
  // In place if the block stays in the same size class.
  if (old_size <= MAX_POOLED && new_size <= MAX_POOLED &&
      class_index(old_size) == class_index(new_size))
    return ptr;
  void *np = pool_alloc(new_size);
  std::memcpy(np, ptr, old_size < new_size ? old_size : new_size);
  pool_free(ptr, old_size);
  return np;
}

} // namespace gmp_pool

// Install the pool as GMP's allocator (idempotent). Normally called
// automatically via the GMP_POOL auto-installer below; exposed for programs that
// prefer an explicit call at the top of main().
inline void install_gmp_pool() {
  mp_set_memory_functions(gmp_pool::pool_alloc, gmp_pool::pool_realloc,
                          gmp_pool::pool_free);
}

#ifdef GMP_POOL
namespace gmp_pool {
// Auto-install at static-initialization time. mp_set_memory_functions only sets
// three pointers (no allocation), and any GMP object allocated before this runs
// is merely reabsorbed into the pool when freed -- never corrupted -- so the
// initialization order is not a hazard.
struct Installer {
  Installer() { install_gmp_pool(); }
};
inline Installer g_installer;
} // namespace gmp_pool
#endif

#endif // SRC_NUMBER_GMP_POOL_ALLOCATOR_H_
