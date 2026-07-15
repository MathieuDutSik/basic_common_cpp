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
//   Two things, both keyed off the single -DGMP_POOL build flag:
//     1. Build with -DGMP_POOL (the Makefiles / CMake set it by default).
//     2. Call one of the installers as the first line of main() (see below). Each
//        expands to a real install under -DGMP_POOL and to nothing otherwise, so
//        the call site stays a single unconditional line and the flag alone turns
//        the pool on/off. Without the call, or without -DGMP_POOL, this is inert.
//
// SERIAL vs THREADED -- CHOOSE THE INSTALLER, NOT A BUILD FLAG
//   The pool exists in two flavours, both always compiled in; which one is active
//   is decided at runtime by which installer main() calls. There is no
//   thread-safety build flag.
//
//   maybe_install_gmp_pool()      -- process-GLOBAL pool state. Fastest (no
//       thread-local indirection, no locks). Use it in single-threaded programs.
//
//   maybe_install_tls_gmp_pool()  -- THREAD-LOCAL pool state: each thread gets its
//       own free-lists and arena. Use it whenever more than one thread of the
//       process may touch GMP concurrently. The global pool is lock-free, so it
//       would corrupt its free-lists under concurrent use; the thread-local pool
//       cannot. Cross-thread frees are safe: a block freed on another thread joins
//       that thread's free-list and is reused there; the memory stays valid.
//
//   Rule of thumb: MPI programs call maybe_install_tls_gmp_pool(). Even though MPI
//   is multi-process (each rank could in principle use the global pool), MPI code
//   -- ours (the dual-description background communication thread) and the MPI
//   runtime itself -- may spawn helper threads that allocate GMP, so the
//   thread-local pool is the safe, take-no-risk default there. Everything else
//   uses maybe_install_gmp_pool().
//
// MEMORY
//   The pool never hands chunks back to the OS during the run -- it grows to the
//   peak working set (the thread-local flavour, to the peak of each live thread).
//   At process exit the OS reclaims everything, so there is no persistent leak (a
//   leak checker may report the pooled chunks as still reachable at exit; that is
//   expected and harmless).
#ifndef SRC_NUMBER_GMP_POOL_ALLOCATOR_H_
#define SRC_NUMBER_GMP_POOL_ALLOCATOR_H_

#include <cstdlib>
#include <cstring>
#include <gmp.h>

namespace gmp_pool {

inline constexpr size_t ALIGN = 16;          // >= alignof(mp_limb_t), >= sizeof(void*)
inline constexpr size_t MAX_POOLED = 4096;   // larger requests go straight to malloc
inline constexpr size_t NCLASS = MAX_POOLED / ALIGN;
inline constexpr size_t CHUNK = size_t(8) << 20; // 8 MB arena chunks

struct FreeNode {
  FreeNode *next;
};

// All the pool's mutable state. Constant-initialized (every member has a constant
// default initializer), so a namespace-scope instance -- global or thread_local --
// needs no dynamic-init guard: the global path stays guard-free and the
// thread_local path pays only the unavoidable TLS access.
struct PoolState {
  FreeNode *freelist[NCLASS] = {};
  char *arena = nullptr;
  size_t arena_left = 0;
};

// 0-based size class for a non-zero size in (0, MAX_POOLED].
inline size_t class_index(size_t size) { return (size + ALIGN - 1) / ALIGN - 1; }

// Core alloc/free/realloc, written once against an explicit state reference; the
// global and thread-local flavours below are thin trampolines over these.
inline void *pool_alloc(PoolState &st, size_t size) {
  if (size == 0)
    size = 1;
  if (size > MAX_POOLED)
    return std::malloc(size);
  size_t ci = class_index(size);
  if (st.freelist[ci] != nullptr) {
    FreeNode *n = st.freelist[ci];
    st.freelist[ci] = n->next;
    return n;
  }
  size_t bytes = (ci + 1) * ALIGN;
  if (st.arena_left < bytes) {
    st.arena = static_cast<char *>(std::malloc(CHUNK));
    st.arena_left = CHUNK;
  }
  void *p = st.arena;
  st.arena += bytes;
  st.arena_left -= bytes;
  return p;
}

inline void pool_free(PoolState &st, void *ptr, size_t size) {
  if (size == 0)
    size = 1;
  if (size > MAX_POOLED) {
    std::free(ptr);
    return;
  }
  size_t ci = class_index(size);
  FreeNode *n = static_cast<FreeNode *>(ptr);
  n->next = st.freelist[ci];
  st.freelist[ci] = n;
}

inline void *pool_realloc(PoolState &st, void *ptr, size_t old_size,
                          size_t new_size) {
  if (old_size == 0)
    old_size = 1;
  if (new_size == 0)
    new_size = 1;
  // In place if the block stays in the same size class.
  if (old_size <= MAX_POOLED && new_size <= MAX_POOLED &&
      class_index(old_size) == class_index(new_size))
    return ptr;
  void *np = pool_alloc(st, new_size);
  std::memcpy(np, ptr, old_size < new_size ? old_size : new_size);
  pool_free(st, ptr, old_size);
  return np;
}

// --- process-global flavour (no locks, no TLS: fastest, single-threaded only) ---
inline PoolState g_global_state;
inline void *global_alloc(size_t size) { return pool_alloc(g_global_state, size); }
inline void *global_realloc(void *ptr, size_t os, size_t ns) {
  return pool_realloc(g_global_state, ptr, os, ns);
}
inline void global_free(void *ptr, size_t size) {
  pool_free(g_global_state, ptr, size);
}

// --- thread-local flavour (one state per thread: safe under concurrent GMP use) --
inline thread_local PoolState g_tls_state;
inline void *tls_alloc(size_t size) { return pool_alloc(g_tls_state, size); }
inline void *tls_realloc(void *ptr, size_t os, size_t ns) {
  return pool_realloc(g_tls_state, ptr, os, ns);
}
inline void tls_free(void *ptr, size_t size) { pool_free(g_tls_state, ptr, size); }

} // namespace gmp_pool

// Install the process-global pool as GMP's allocator (idempotent). Call once, at
// the top of main(), before any GMP object is created. Single-threaded programs.
inline void install_gmp_pool() {
  mp_set_memory_functions(gmp_pool::global_alloc, gmp_pool::global_realloc,
                          gmp_pool::global_free);
}

// Install the thread-local pool as GMP's allocator (idempotent). Call once, at the
// top of main(). Use in any program where several threads may use GMP at once
// (MPI programs -- take-no-risk default).
inline void install_tls_gmp_pool() {
  mp_set_memory_functions(gmp_pool::tls_alloc, gmp_pool::tls_realloc,
                          gmp_pool::tls_free);
}

// The intended entry points: put one of these as the first line of main(). Each
// installs the corresponding pool when built with -DGMP_POOL and is a no-op
// otherwise, so the single GMP_POOL build flag turns the pool on/off with no
// #ifdef at the call site. Use the plain one for serial programs and the _tls_ one
// for MPI / multi-threaded programs.
inline void maybe_install_gmp_pool() {
#ifdef GMP_POOL
  install_gmp_pool();
#endif
}

inline void maybe_install_tls_gmp_pool() {
#ifdef GMP_POOL
  install_tls_gmp_pool();
#endif
}

#endif // SRC_NUMBER_GMP_POOL_ALLOCATOR_H_
