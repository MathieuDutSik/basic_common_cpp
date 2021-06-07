#ifndef INCLUDE_GRAPH_TRACES
#define INCLUDE_GRAPH_TRACES

#include "traces.h"
#include <vector>
#include <iostream>
#include <ExceptionEnding.h>

struct DataTraces {
public:
  int n;
  int nbAdjacent;
  int* lab1;
  int* ptn;
  int* orbits;
  sparsegraph sg1;
  sparsegraph cg1;
  DataTraces(int _n, int _nbAdjacent) : n(_n), nbAdjacent(_nbAdjacent)
  {
    lab1 = (int*)malloc(n * sizeof(int));
    ptn = (int*)malloc(n * sizeof(int));
    orbits = (int*)malloc(n * sizeof(int));
    orbits = (int*)malloc(n * sizeof(int));
    //
    sg1.nv = n;
    sg1.nde = nbAdjacent;
    sg1.v = (size_t*)malloc(n * sizeof(size_t));
    sg1.d = (int*)malloc(n * sizeof(int));
    sg1.e = (int*)malloc(nbAdjacent * sizeof(int));
    sg1.w = NULL;
    sg1.vlen = n;
    sg1.dlen = n;
    sg1.elen = nbAdjacent;
    sg1.wlen = 0;
    //
    cg1.d = NULL; // set to NULL, but it may or may not be expanded.
    cg1.v = NULL;
    cg1.e = NULL;
    cg1.w = NULL;
  }
  ~DataTraces()
  {
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
};


template<typename Tidx>
std::vector<Tidx> TRACES_GetCanonicalOrdering_Arr(DataTraces* DT)
{
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    int n = DT->n;

    options.getcanon = TRUE;
    options.defaultptn = FALSE;

    Traces(&DT->sg1, DT->lab1, DT->ptn, DT->orbits, &options, &stats, &DT->cg1);
    std::vector<Tidx> V(n);
    for (int i=0; i<n; i++)
      V[DT->lab1[i]] = i;
    return V;
}


template<typename Tgr>
DataTraces GetDataTraces_from_G(Tgr const& eGR)
{
  int n = eGR.GetNbVert();
  int nbAdjacent = eGR.GetNbAdjacent();
  bool HasVertexColor = eGR.GetHasVertexColor();
  if (!HasVertexColor) {
    std::cerr << "The graph should have vertex color\n";
    throw TerminalException{1};
  }
  // allocating the DataTraces
  DataTraces DT(n, nbAdjacent);
  // Determining the color symbolic information
  int numcells=0;
  for (int i=0; i<n; i++) {
    int eVal = 1 + eGR.GetColor(i);
    if (eVal > numcells)
      numcells = eVal;
  }
  std::vector<int> ListPartSize(numcells,0);
  for (int i=0; i<n; i++)
    ListPartSize[eGR.GetColor(i)]++;
  std::vector<int> ListShift(numcells,0);
  for (int icell=1; icell<numcells; icell++)
    ListShift[icell] = ListShift[icell-1] + ListPartSize[icell-1];
  // lab1 construction
  for (int i=0; i<n; i++) {
    int icell = eGR.GetColor(i);
    DT.lab1[ListShift[icell]] = i;
    ListShift[icell]++;
  }
  // ptn construction
  for (int i=0; i<n; i++) DT.ptn[i] = NAUTY_INFINITY;
  for (int icell=0; icell<numcells; icell++)
    DT.ptn[ListShift[icell] - 1] = 0;
  // The adjacencies
  int pos = 0;
  for (int i=0; i<n; i++) {
    std::vector<int> LAdj = eGR.Adjacency(i);
    int len = LAdj.size();
    DT.sg1.d[i] = len;
    DT.sg1.v[i] = pos;
    for (auto & eAdj : LAdj) {
      DT.sg1.e[pos] = eAdj;
      pos++;
    }
  }
  // Returning the final entry
  return DT;
}





template<typename Tgr, typename Tidx>
std::vector<Tidx> TRACES_GetCanonicalOrdering(Tgr const& eGR)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    /* Declare and initialize sparse graph structures */
    SG_DECL(sg1);
    SG_DECL(cg1);

    /* Reading key graph variables */
    int n = eGR.GetNbVert();
    int nbAdjacent = eGR.GetNbAdjacent();
    bool HasVertexColor = eGR.GetHasVertexColor();

    /* Select option for canonical labelling */
    options.getcanon = TRUE;

    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    if (HasVertexColor) {
      options.defaultptn = FALSE;
      int numcells=0;
      for (int i=0; i<n; i++) {
        int eVal = 1 + eGR.GetColor(i);
        if (eVal > numcells)
          numcells = eVal;
      }
      std::vector<int> ListPartSize(numcells,0);
      for (int i=0; i<n; i++)
        ListPartSize[eGR.GetColor(i)]++;
      std::vector<int> ListShift(numcells,0);
      for (int icell=1; icell<numcells; icell++)
        ListShift[icell] = ListShift[icell-1] + ListPartSize[icell-1];
      // lab1 construction
      for (int i=0; i<n; i++) {
        int icell = eGR.GetColor(i);
        lab1[ListShift[icell]] = i;
        ListShift[icell]++;
      }
      // ptn construction
      for (int i=0; i<n; i++) ptn[i] = NAUTY_INFINITY;
      for (int icell=0; icell<numcells; icell++)
        ptn[ListShift[icell] - 1] = 0;
    }

    /* Now make the graph */
    SG_ALLOC(sg1,n,nbAdjacent,"malloc");
    sg1.nv = n;              /* Number of vertices */
    sg1.nde = nbAdjacent;           /* Number of directed edges */

    int pos = 0;
    for (int i=0; i<n; i++) {
      std::vector<int> LAdj = eGR.Adjacency(i);
      int len = LAdj.size();
      sg1.d[i] = len;
      sg1.v[i] = pos;
      for (auto & eAdj : LAdj) {
        sg1.e[pos] = eAdj;
        pos++;
      }
    }

    Traces(&sg1,lab1,ptn,orbits,&options,&stats,&cg1);
    std::vector<Tidx> V(n);
    for (int i=0; i<n; i++)
      V[lab1[i]] = i;

    DYNFREE(lab1,lab1_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    SG_FREE(sg1);
    SG_FREE(cg1);
    return V;
}



template<typename Tgr, typename Tidx>
std::vector<Tidx> TRACES_GetCanonicalOrdering_Arr_Test(Tgr const& eGR)
{
  DataTraces DT = GetDataTraces_from_G(eGR);
  std::cerr << "After TRACES_GetCanonicalOrdering_Arr_Test\n";
  return TRACES_GetCanonicalOrdering_Arr<Tidx>(&DT);
}




template<typename Tidx>
std::vector<std::vector<Tidx>> TRACES_GetListGenerators_Arr(DataTraces* DT, int n_last)
{
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    permnode *gens;
    options.generators = &gens;
    gens = NULL;
    freeschreier(NULL,&gens);
    options.getcanon = FALSE;

    /*
    int n = DT->n;
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    */

    options.defaultptn = FALSE;

    Traces(&DT->sg1, DT->lab1, DT->ptn, DT->orbits, &options, &stats, NULL);

    std::vector<std::vector<Tidx>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<Tidx> V(n_last);
          for (int i=0; i<n_last; i++)
            V[i] = pn->p[i];
          ListGen.push_back(V);
          //
          pn = pn->next;
        } while (pn != gens);
    }
    freeschreier(NULL,&gens);
    schreier_freedyn();
    return ListGen;
}


template<typename Tgr, typename Tidx>
std::vector<std::vector<Tidx>> TRACES_GetListGenerators(Tgr const& eGR, int n_last)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    /* Declare the generator stuff */
    permnode *gens;
    options.generators = &gens;
    gens = NULL;
    freeschreier(NULL,&gens);
    options.getcanon = FALSE;

    /* Reading key graph variables */
    int n = eGR.GetNbVert();
    int nbAdjacent = eGR.GetNbAdjacent();
    bool HasVertexColor = eGR.GetHasVertexColor();

    /* Declare and initialize sparse graph structures */
    SG_DECL(sg1);
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    /*
      "lab" and "ptn" contain the information of the initial partition.
      It is a complex construction that goes deep into the internals of nauty/traces.
      ---Initially ptn contains only two possible values: NAUTY_INFINITY and 0
     */
    if (HasVertexColor) {
      options.defaultptn = FALSE;
      int numcells=0;
      for (int i=0; i<n; i++) {
        int eVal = 1 + eGR.GetColor(i);
        if (eVal > numcells)
          numcells = eVal;
      }
      std::vector<int> ListPartSize(numcells,0);
      for (int i=0; i<n; i++)
        ListPartSize[eGR.GetColor(i)]++;
      std::vector<int> ListShift(numcells,0);
      for (int icell=1; icell<numcells; icell++)
        ListShift[icell] = ListShift[icell-1] + ListPartSize[icell-1];
      // lab1 construction
      for (int i=0; i<n; i++) {
        int icell = eGR.GetColor(i);
        lab1[ListShift[icell]] = i;
        ListShift[icell]++;
      }
      // ptn construction
      for (int i=0; i<n; i++) ptn[i] = NAUTY_INFINITY;
      for (int icell=0; icell<numcells; icell++)
        ptn[ListShift[icell] - 1] = 0;
    }

    /* Now make the graph */
    SG_ALLOC(sg1,n,nbAdjacent,"malloc");
    sg1.nv = n;              /* Number of vertices */
    sg1.nde = nbAdjacent;           /* Number of directed edges */
    int pos = 0;
    for (int i=0; i<n; i++) {
      std::vector<int> LAdj = eGR.Adjacency(i);
      int len = LAdj.size();
      sg1.d[i] = len;
      sg1.v[i] = pos;
      for (auto & eAdj : LAdj) {
        sg1.e[pos] = eAdj;
        pos++;
      }
    }
    /* Calling Traces */
    Traces(&sg1,lab1,ptn,orbits,&options,&stats,NULL);
    /* Extracting the list of generators */
    std::vector<std::vector<Tidx>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<Tidx> V(n_last);
          for (int i=0; i<n_last; i++)
            V[i] = pn->p[i];
          ListGen.push_back(V);
          //
          pn = pn->next;
        } while (pn != gens);
    }
    freeschreier(NULL,&gens);
    schreier_freedyn();

    DYNFREE(lab1,lab1_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    SG_FREE(sg1);
    return ListGen;
}


template<typename Tgr, typename Tidx>
std::vector<std::vector<Tidx>> TRACES_GetListGenerators_Arr_Test(Tgr const& eGR, int n_last)
{
  DataTraces DT = GetDataTraces_from_G(eGR);
  std::cerr << "After TRACES_GetListGenerators_Arr_Test\n";
  return TRACES_GetListGenerators_Arr<Tidx>(&DT, n_last);
}


template<typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> TRACES_GetCanonicalOrdering_ListGenerators_Arr(DataTraces* DT, int n_last)
{
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    permnode *gens;
    options.generators = &gens;
    gens = NULL;
    freeschreier(NULL,&gens);

    int n = DT->n;
    options.getcanon = TRUE;
    options.defaultptn = FALSE;

    Traces(&DT->sg1, DT->lab1, DT->ptn, DT->orbits, &options, &stats, &DT->cg1);
    std::vector<Tidx> V(n);
    for (int i=0; i<n; i++)
      V[DT->lab1[i]] = i;
    std::vector<std::vector<Tidx>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<Tidx> V(n_last);
          for (int i=0; i<n_last; i++)
            V[i] = pn->p[i];
          ListGen.push_back(V);
          //
          pn = pn->next;
        } while (pn != gens);
    }
    freeschreier(NULL,&gens);
    schreier_freedyn();
    return {std::move(V), std::move(ListGen)};
}



template<typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> TRACES_GetCanonicalOrdering_ListGenerators(Tgr const& eGR, int n_last)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_TRACES(options);
    TracesStats stats;
    /* Declare the generator stuff */
    permnode *gens;
    options.generators = &gens;
    gens = NULL;
    freeschreier(NULL,&gens);

    /* Declare and initialize sparse graph structures */
    SG_DECL(sg1);
    SG_DECL(cg1);

    /* Reading key graph variables */
    int n = eGR.GetNbVert();
    int nbAdjacent = eGR.GetNbAdjacent();
    bool HasVertexColor = eGR.GetHasVertexColor();

    /* Select option for canonical labelling */
    options.getcanon = TRUE;

    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    if (HasVertexColor) {
      options.defaultptn = FALSE;
      int numcells=0;
      for (int i=0; i<n; i++) {
        int eVal = 1 + eGR.GetColor(i);
        if (eVal > numcells)
          numcells = eVal;
      }
      std::vector<int> ListPartSize(numcells,0);
      for (int i=0; i<n; i++)
        ListPartSize[eGR.GetColor(i)]++;
      std::vector<int> ListShift(numcells,0);
      for (int icell=1; icell<numcells; icell++)
        ListShift[icell] = ListShift[icell-1] + ListPartSize[icell-1];
      // lab1 construction
      for (int i=0; i<n; i++) {
        int icell = eGR.GetColor(i);
        lab1[ListShift[icell]] = i;
        ListShift[icell]++;
      }
      // ptn construction
      for (int i=0; i<n; i++) ptn[i] = NAUTY_INFINITY;
      for (int icell=0; icell<numcells; icell++)
        ptn[ListShift[icell] - 1] = 0;
    }

    /* Now make the graph */
    SG_ALLOC(sg1,n,nbAdjacent,"malloc");
    sg1.nv = n;              /* Number of vertices */
    sg1.nde = nbAdjacent;           /* Number of directed edges */

    int pos = 0;
    for (int i=0; i<n; i++) {
      std::vector<int> LAdj = eGR.Adjacency(i);
      int len = LAdj.size();
      sg1.d[i] = len;
      sg1.v[i] = pos;
      for (auto & eAdj : LAdj) {
        sg1.e[pos] = eAdj;
        pos++;
      }
    }

    Traces(&sg1,lab1,ptn,orbits,&options,&stats,&cg1);
    // Extracting the canonical ordering
    std::vector<Tidx> V(n);
    for (int i=0; i<n; i++)
      V[lab1[i]] = i;
    // Extracting the list of generators
    std::vector<std::vector<Tidx>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<Tidx> V(n_last);
          for (int i=0; i<n_last; i++)
            V[i] = pn->p[i];
          ListGen.push_back(V);
          //
          pn = pn->next;
        } while (pn != gens);
    }
    freeschreier(NULL,&gens);
    schreier_freedyn();

    DYNFREE(lab1,lab1_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    SG_FREE(sg1);
    SG_FREE(cg1);
    return {std::move(V), std::move(ListGen)};
}


#endif
