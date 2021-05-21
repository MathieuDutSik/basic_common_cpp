#ifndef INCLUDE_GRAPH_TRACES
#define INCLUDE_GRAPH_TRACES

#include "traces.h"
#include <vector>
#include <iostream>

template<typename Tgr>
std::vector<unsigned int> TRACES_GetCanonicalOrdering(Tgr const& eGR)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
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
    DYNALLOC1(int,map,map_sz,n,"malloc");
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
    std::vector<unsigned int> V(n);
    for (int i=0; i<n; i++)
      V[lab1[i]] = i;

    DYNFREE(lab1,lab1_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    DYNFREE(map,map_sz);
    SG_FREE(sg1);
    SG_FREE(cg1);
    return V;
}


template<typename Tgr>
std::vector<std::vector<unsigned int>> TRACES_GetListGenerators(Tgr const& eGR)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
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
    DYNALLOC1(int,map,map_sz,n,"malloc");
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
    std::vector<std::vector<unsigned int>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<unsigned int> V(n);
          for (int i=0; i<n; i++)
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
    DYNFREE(map,map_sz);
    SG_FREE(sg1);
    return ListGen;
}



template<typename Tgr>
std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> TRACES_GetCanonicalOrdering_ListGenerators(Tgr const& eGR)
{
    DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
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
    DYNALLOC1(int,map,map_sz,n,"malloc");
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
    std::vector<unsigned int> V(n);
    for (int i=0; i<n; i++)
      V[lab1[i]] = i;
    // Extracting the list of generators
    std::vector<std::vector<unsigned int>> ListGen;
    if (gens) {
      permnode* pn = gens;
      do
        {
          std::vector<unsigned int> V(n);
          for (int i=0; i<n; i++)
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
    DYNFREE(map,map_sz);
    SG_FREE(sg1);
    SG_FREE(cg1);
    return {std::move(V), std::move(ListGen)};
}


#endif
