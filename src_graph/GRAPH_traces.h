// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_TRACES_H_
#define SRC_GRAPH_GRAPH_TRACES_H_

#include "ExceptionsFunc.h"
#include "Timings.h"
#include "traces.h"
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#ifdef TIMINGS
#define TIMINGS_TRACES
#endif

struct DataTraces {
public:
  size_t n;
  size_t nbAdjacent;
  int *lab1;
  int *ptn;
  int *orbits;
  sparsegraph sg1;
  sparsegraph cg1;
  DataTraces(size_t _n, size_t _nbAdjacent) : n(_n), nbAdjacent(_nbAdjacent) {
    lab1 = reinterpret_cast<int *>(malloc(n * sizeof(int)));
    ptn = reinterpret_cast<int *>(malloc(n * sizeof(int)));
    orbits = reinterpret_cast<int *>(malloc(n * sizeof(int)));
    // We allocate the arrays for Traces
    sg1.nv = static_cast<int>(n);
    sg1.nde = static_cast<int>(nbAdjacent);
    sg1.v = reinterpret_cast<size_t *>(malloc(n * sizeof(size_t)));
    sg1.d = reinterpret_cast<int *>(malloc(n * sizeof(int)));
    sg1.e = reinterpret_cast<int *>(malloc(nbAdjacent * sizeof(int)));
    sg1.w = NULL;
    sg1.vlen = static_cast<int>(n);
    sg1.dlen = static_cast<int>(n);
    sg1.elen = static_cast<int>(nbAdjacent);
    sg1.wlen = 0;
    // We imitate SG_DECL here
    cg1.nv = 0;
    cg1.nde = 0;
    // set to NULL, but it may or may not be expanded.
    cg1.d = NULL;
    cg1.v = NULL;
    cg1.e = NULL;
    cg1.w = NULL;
    cg1.vlen = 0;
    cg1.dlen = 0;
    cg1.elen = 0;
    cg1.wlen = 0;
  }
  ~DataTraces() {
    free(lab1);
    free(ptn);
    free(orbits);
    //
    free(sg1.d);
    free(sg1.v);
    free(sg1.e);
    //
    free(cg1.d);
    free(cg1.v);
    free(cg1.e);
  }
  DataTraces() = delete;
  DataTraces(const DataTraces &) = delete;
  DataTraces(DataTraces &&x) {
    n = x.n;
    nbAdjacent = x.nbAdjacent;
    //
    lab1 = x.lab1;
    ptn = x.ptn;
    orbits = x.orbits;
    x.lab1 = NULL;
    x.ptn = NULL;
    x.orbits = NULL;
    //
    sg1.nv = x.sg1.nv;
    sg1.nde = x.sg1.nde;
    sg1.vlen = x.sg1.vlen;
    sg1.dlen = x.sg1.dlen;
    sg1.elen = x.sg1.elen;
    sg1.d = x.sg1.d;
    sg1.v = x.sg1.v;
    sg1.e = x.sg1.e;
    x.sg1.d = NULL;
    x.sg1.v = NULL;
    x.sg1.e = NULL;
    //
    //
    cg1.nv = x.cg1.nv;
    cg1.nde = x.cg1.nde;
    cg1.vlen = x.cg1.vlen;
    cg1.dlen = x.cg1.dlen;
    cg1.elen = x.cg1.elen;
    cg1.d = x.cg1.d;
    cg1.v = x.cg1.v;
    cg1.e = x.cg1.e;
    x.cg1.d = NULL;
    x.cg1.v = NULL;
    x.cg1.e = NULL;
  }
  DataTraces &operator=(const DataTraces &) = delete;
};

template <typename Tidx> void TRACES_LimitCheck(size_t n) {
  if (n >= size_t(std::numeric_limits<Tidx>::max())) {
    std::cerr << "Error in TRACES_LimitCheck\n";
    std::cerr << "We have n=" << n << " std::numeric_limits<Tidx>::max()="
              << std::numeric_limits<Tidx>::max() << "\n";
    throw TerminalException{1};
  }
}

template <typename Tidx>
std::vector<Tidx>
TRACES_GetCanonicalOrdering_Arr(DataTraces &DT,
                                [[maybe_unused]] std::ostream &os) {
  size_t n = DT.n;
  TRACES_LimitCheck<Tidx>(n);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;

  options.getcanon = TRUE;
  options.defaultptn = FALSE;

  Traces(&DT.sg1, DT.lab1, DT.ptn, DT.orbits, &options, &stats, &DT.cg1);
  std::vector<Tidx> V(n);
  for (size_t i = 0; i < n; i++)
    V[DT.lab1[i]] = Tidx(i);
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetCanonicalOrdering_Arr|=" << time << "\n";
#endif
  return V;
}

template <typename Tgr> void Assign_sg(Tgr const &eGR, sparsegraph *sg) {
  size_t n = eGR.GetNbVert();
  size_t pos = 0;
  for (size_t i = 0; i < n; i++) {
    std::vector<size_t> LAdj = eGR.Adjacency(i);
    size_t len = LAdj.size();
    sg->d[i] = static_cast<int>(len);
    sg->v[i] = static_cast<int>(pos);
    for (auto &eAdj : LAdj) {
      sg->e[pos] = static_cast<int>(eAdj);
      pos++;
    }
  }
}

template <typename Tgr>
void Assign_lab1_ptn(Tgr const &eGR, int *lab1, int *ptn) {
  size_t n = eGR.GetNbVert();
  size_t numcells = 0;
  for (size_t i = 0; i < n; i++) {
    size_t eVal = 1 + eGR.GetColor(i);
    if (eVal > numcells)
      numcells = eVal;
  }
  std::vector<int> ListPartSize(numcells, 0);
  for (size_t i = 0; i < n; i++)
    ListPartSize[eGR.GetColor(i)]++;
  std::vector<size_t> ListShift(numcells, 0);
  for (size_t icell = 1; icell < numcells; icell++)
    ListShift[icell] = ListShift[icell - 1] + ListPartSize[icell - 1];
  // lab1 construction
  for (size_t i = 0; i < n; i++) {
    size_t icell = eGR.GetColor(i);
    lab1[ListShift[icell]] = static_cast<int>(i);
    ListShift[icell]++;
  }
  // ptn construction
  for (size_t i = 0; i < n; i++)
    ptn[i] = NAUTY_INFINITY;
  for (size_t icell = 0; icell < numcells; icell++)
    ptn[ListShift[icell] - 1] = 0;
}

template <typename Tgr>
void Assign_lab1_ptn_sg(Tgr const &eGR, bool const &HasVertexColor, int *lab1,
                        int *ptn, sparsegraph *sg) {
  if (HasVertexColor) {
    Assign_lab1_ptn(eGR, lab1, ptn);
  }
  Assign_sg(eGR, sg);
}

template <typename Tgr> DataTraces *GetDataTraces_from_G(Tgr const &eGR) {
  size_t n = eGR.GetNbVert();
  size_t nbAdjacent = eGR.GetNbAdjacent();
  bool HasVertexColor = eGR.GetHasVertexColor();
  if (!HasVertexColor) {
    std::cerr << "The graph should have vertex color\n";
    throw TerminalException{1};
  }
  // allocating the DataTraces
  DataTraces *DT = new DataTraces(n, nbAdjacent);
  // Determining the color symbolic information
  Assign_lab1_ptn_sg(eGR, HasVertexColor, DT->lab1, DT->ptn, &(DT->sg1));
  // Returning the final entry
  return DT;
}

template <typename Tgr, typename Tidx>
std::vector<Tidx>
TRACES_GetCanonicalOrdering(Tgr const &eGR, [[maybe_unused]] std::ostream &os) {
  size_t n = eGR.GetNbVert();
  TRACES_LimitCheck<Tidx>(n);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  DYNALLSTAT(int, lab1, lab1_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  /* Declare and initialize sparse graph structures */
  SG_DECL(sg1);
  SG_DECL(cg1);

  /* Reading key graph variables */
  size_t nbAdjacent = eGR.GetNbAdjacent();
  bool HasVertexColor = eGR.GetHasVertexColor();

  /* Now make the graph */
  int n_i = static_cast<int>(n);
  int nbAdjacent_i = static_cast<int>(nbAdjacent);
  SG_ALLOC(sg1, n_i, nbAdjacent_i, "malloc");
  sg1.nv = n_i;           /* Number of vertices */
  sg1.nde = nbAdjacent_i; /* Number of directed edges */

  /* Select option for canonical labelling */
  options.getcanon = TRUE;

  int m = SETWORDSNEEDED(n_i);
  nauty_check(WORDSIZE, m, n_i, NAUTYVERSIONID);

  DYNALLOC1(int, lab1, lab1_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

  if (HasVertexColor)
    options.defaultptn = FALSE;
  Assign_lab1_ptn_sg(eGR, HasVertexColor, lab1, ptn, &sg1);

  Traces(&sg1, lab1, ptn, orbits, &options, &stats, &cg1);
  std::vector<Tidx> V(n);
  for (size_t i = 0; i < n; i++)
    V[lab1[i]] = Tidx(i);

  DYNFREE(lab1, lab1_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  SG_FREE(sg1);
  SG_FREE(cg1);
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetCanonicalOrdering|=" << time << "\n";
#endif
  return V;
}

template <typename Tgr, typename Tidx>
std::vector<Tidx> TRACES_GetCanonicalOrdering_Arr_Test(Tgr const &eGR) {
  DataTraces *DT = GetDataTraces_from_G(eGR);
  std::vector<Tidx> V = TRACES_GetCanonicalOrdering_Arr<Tidx>(*DT, std::cerr);
  delete DT;
  return V;
}

template <typename Tidx>
std::vector<std::vector<Tidx>>
TRACES_GetListGenerators_Arr(DataTraces &DT, size_t const &n_last,
                             [[maybe_unused]] std::ostream &os) {
  TRACES_LimitCheck<Tidx>(n_last);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  permnode *gens;
  options.generators = &gens;
  gens = NULL;
  freeschreier(NULL, &gens);
  options.getcanon = FALSE;

  options.defaultptn = FALSE;

  Traces(&DT.sg1, DT.lab1, DT.ptn, DT.orbits, &options, &stats, NULL);

  std::vector<std::vector<Tidx>> ListGen;
  if (gens) {
    permnode *pn = gens;
    std::vector<Tidx> V(n_last);
    do {
      for (size_t i = 0; i < n_last; i++)
        V[i] = pn->p[i];
      ListGen.push_back(V);
      //
      pn = pn->next;
    } while (pn != gens);
  }
  freeschreier(NULL, &gens);
  schreier_freedyn();
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetListGenerators_Arr|=" << time << "\n";
#endif
  return ListGen;
}

template <typename Tidx>
void ReadListGen(permnode *gens, std::vector<std::vector<Tidx>> &ListGen,
                 size_t const &n_last) {
  if (gens) {
    permnode *pn = gens;
    do {
      std::vector<Tidx> V(n_last);
      for (size_t i = 0; i < n_last; i++)
        V[i] = pn->p[i];
      ListGen.push_back(V);
      //
      pn = pn->next;
    } while (pn != gens);
  }
}

template <typename Tgr, typename Tidx>
std::vector<std::vector<Tidx>>
TRACES_GetListGenerators(Tgr const &eGR, size_t const &n_last,
                         [[maybe_unused]] std::ostream &os) {
  TRACES_LimitCheck<Tidx>(n_last);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  DYNALLSTAT(int, lab1, lab1_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  /* Declare the generator stuff */
  permnode *gens;
  options.generators = &gens;
  gens = NULL;
  freeschreier(NULL, &gens);
  options.getcanon = FALSE;

  /* Reading key graph variables */
  size_t n = eGR.GetNbVert();
  size_t nbAdjacent = eGR.GetNbAdjacent();
  bool HasVertexColor = eGR.GetHasVertexColor();

  /* Declare and initialize sparse graph structures */
  SG_DECL(sg1);
  int n_i = static_cast<int>(n);
  int m = SETWORDSNEEDED(n_i);
  nauty_check(WORDSIZE, m, n_i, NAUTYVERSIONID);

  DYNALLOC1(int, lab1, lab1_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
  /*
    "lab" and "ptn" contain the information of the initial partition.
    It is a complex construction that goes deep into the internals of
    nauty/traces.
    ---Initially ptn contains only two possible values: NAUTY_INFINITY and 0
  */
  if (HasVertexColor) {
    options.defaultptn = FALSE;
    Assign_lab1_ptn(eGR, lab1, ptn);
  }

  /* Now make the graph */
  int nbAdjacent_i = static_cast<int>(nbAdjacent);
  SG_ALLOC(sg1, n, nbAdjacent, "malloc");
  sg1.nv = n_i;           /* Number of vertices */
  sg1.nde = nbAdjacent_i; /* Number of directed edges */
  Assign_sg(eGR, &sg1);
  /* Calling Traces */
  Traces(&sg1, lab1, ptn, orbits, &options, &stats, NULL);
  /* Extracting the list of generators */
  std::vector<std::vector<Tidx>> ListGen;
  ReadListGen(gens, ListGen, n_last);
  freeschreier(NULL, &gens);
  schreier_freedyn();

  DYNFREE(lab1, lab1_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  SG_FREE(sg1);
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetListGenerators|=" << time << "\n";
#endif
  return ListGen;
}

template <typename Tgr, typename Tidx>
std::vector<std::vector<Tidx>>
TRACES_GetListGenerators_Arr_Test(Tgr const &eGR, size_t const &n_last) {
  DataTraces *DT = GetDataTraces_from_G(eGR);
  std::vector<std::vector<Tidx>> ret =
      TRACES_GetListGenerators_Arr<Tidx>(*DT, n_last, std::cerr);
  delete DT;
  return ret;
}

template <typename TidxC, typename TidxG>
std::pair<std::vector<TidxC>, std::vector<std::vector<TidxG>>>
TRACES_GetCanonicalOrdering_ListGenerators_Arr(
    DataTraces &DT, size_t const &n_last, [[maybe_unused]] std::ostream &os) {
  size_t n = size_t(DT.n);
  TRACES_LimitCheck<TidxC>(n);
  TRACES_LimitCheck<TidxG>(n_last);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  permnode *gens;
  options.generators = &gens;
  gens = NULL;
  freeschreier(NULL, &gens);

  options.getcanon = TRUE;
  options.defaultptn = FALSE;

  Traces(&DT.sg1, DT.lab1, DT.ptn, DT.orbits, &options, &stats, &DT.cg1);
  std::vector<TidxC> V(n);
  for (size_t i = 0; i < n; i++)
    V[DT.lab1[i]] = TidxC(i);
  std::vector<std::vector<TidxG>> ListGen;
  ReadListGen(gens, ListGen, n_last);
  freeschreier(NULL, &gens);
  schreier_freedyn();
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetCanonicalOrdering_ListGenerators_Arr|=" << time << "\n";
#endif
  return {std::move(V), std::move(ListGen)};
}

template <typename Tgr, typename TidxC, typename TidxG>
std::pair<std::vector<TidxC>, std::vector<std::vector<TidxG>>>
TRACES_GetCanonicalOrdering_ListGenerators(Tgr const &eGR, size_t n_last,
                                           [[maybe_unused]] std::ostream &os) {
  size_t n = eGR.GetNbVert();
  TRACES_LimitCheck<TidxC>(n);
  TRACES_LimitCheck<TidxG>(n_last);
#ifdef TIMINGS_TRACES
  MicrosecondTime time;
#endif
  DYNALLSTAT(int, lab1, lab1_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  /* Declare the generator stuff */
  permnode *gens;
  options.generators = &gens;
  gens = NULL;
  freeschreier(NULL, &gens);

  /* Declare and initialize sparse graph structures */
  SG_DECL(sg1);
  SG_DECL(cg1);

  /* Reading key graph variables */
  int nbAdjacent = eGR.GetNbAdjacent();
  bool HasVertexColor = eGR.GetHasVertexColor();

  /* Select option for canonical labelling */
  options.getcanon = TRUE;

  int n_i = static_cast<int>(n);
  int m = SETWORDSNEEDED(n_i);
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  DYNALLOC1(int, lab1, lab1_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
  if (HasVertexColor) {
    options.defaultptn = FALSE;
    int numcells = 0;
    for (size_t i = 0; i < n; i++) {
      int eVal = 1 + eGR.GetColor(i);
      if (eVal > numcells)
        numcells = eVal;
    }
    std::vector<int> ListPartSize(numcells, 0);
    for (size_t i = 0; i < n; i++)
      ListPartSize[eGR.GetColor(i)]++;
    std::vector<int> ListShift(numcells, 0);
    for (int icell = 1; icell < numcells; icell++)
      ListShift[icell] = ListShift[icell - 1] + ListPartSize[icell - 1];
    // lab1 construction
    for (size_t i = 0; i < n; i++) {
      int icell = eGR.GetColor(i);
      lab1[ListShift[icell]] = i;
      ListShift[icell]++;
    }
    // ptn construction
    for (size_t i = 0; i < n; i++)
      ptn[i] = NAUTY_INFINITY;
    for (int icell = 0; icell < numcells; icell++)
      ptn[ListShift[icell] - 1] = 0;
  }

  /* Now make the graph */
  SG_ALLOC(sg1, n, nbAdjacent, "malloc");
  sg1.nv = n;           /* Number of vertices */
  sg1.nde = nbAdjacent; /* Number of directed edges */
  Assign_sg(eGR, &sg1);
  Traces(&sg1, lab1, ptn, orbits, &options, &stats, &cg1);
  // Extracting the canonical ordering
  std::vector<TidxC> V(n);
  for (size_t i = 0; i < n; i++)
    V[lab1[i]] = i;
  // Extracting the list of generators
  std::vector<std::vector<TidxG>> ListGen;
  ReadListGen(gens, ListGen, n_last);
  freeschreier(NULL, &gens);
  schreier_freedyn();

  DYNFREE(lab1, lab1_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  SG_FREE(sg1);
  SG_FREE(cg1);
#ifdef TIMINGS_TRACES
  os << "|TRA: TRACES_GetCanonicalOrdering_ListGenerators|=" << time << "\n";
#endif
  return {std::move(V), std::move(ListGen)};
}

// clang-format off
#endif  // SRC_GRAPH_GRAPH_TRACES_H_
// clang-format on
