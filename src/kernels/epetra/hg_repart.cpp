/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
// this file is an adaptation of an example from Trilinos 11.12.1, isorropia package.
// original header is below.

//@HEADER
//
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//@HEADER

//--------------------------------------------------------------------
//This file is a self-contained example of creating an Epetra_RowMatrix
//object, and using Isorropia to create a rebalanced copy of it using
//Zoltan's Hypergraph partitioning.
//Hypergraph edge weights are used to influence the repartitioning.
//--------------------------------------------------------------------


#include "phist_config.h"

#include "phist_trilinos_macros.h"
#include "phist_tools.h"

#include "Teuchos_StandardCatchMacros.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>


#ifndef PHIST_HAVE_ISORROPIA

void repartition(Teuchos::RCP<const Epetra_CrsMatrix> rowmatrix,
                 Teuchos::RCP<      Epetra_Map>& bal_map,
                 Teuchos::RCP<      Epetra_CrsMatrix>& bal_mat, bool redist,
                 int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}
#else
//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_Exception.hpp>

// tools to analyze the quality of the partitioning
//#include "ispatest_lbeval_utils.hpp"

namespace phist 
{
  namespace epetra_internal 
  {

    // on input, rowmatrix must be Filled(). on output, bal_map points to the new (repartitioned) map,
    // and bal_mat contains the repartitioned matrix *if redist==true*. Otherwise, bal_mat=Teuchos::null
    // on output.
    void repartition(Teuchos::RCP<const Epetra_RowMatrix> rowmatrix,
                     Teuchos::RCP<      Epetra_Map>& bal_map,
                     Teuchos::RCP<      Epetra_CrsMatrix>& bal_matrix, bool redist,
                     int *iflag)
    {

      //We'll need a Teuchos::ParameterList object to pass to the
      //Isorropia::Epetra::Partitioner class.
      Teuchos::ParameterList paramlist;

#ifdef HAVE_ISORROPIA_ZOLTAN
      // If Zoltan is available, we'll specify that the Zoltan package be
      // used for the partitioning operation, by creating a parameter
      // sublist named "Zoltan".
      // In the sublist, we'll set parameters that we want sent to Zoltan.
      // (As it turns out, Isorropia selects Zoltan's hypergraph partitioner
      //  by default, so we don't actually need to specify it. But it's
      //  useful for illustration...)

      paramlist.set("PARTITIONING METHOD", "HYPERGRAPH");

#else
      // If Zoltan is not available, a simple linear partitioner will be
      // used to partition such that the number of nonzeros is equal (or
      // close to equal) on each processor. No parameter is necessary to
      // specify this.
# ifndef I_WARNED_YOU
# define I_WARNED_YOU
# warning "compiling with Isorropia but without Zoltan, repartitioning will fallback to a simple linear partitioner"
# endif
#endif
      Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner=Teuchos::null;
      
      // check for the combination 64-bit gids + Isorropia, this tends to cause an exception
      // even if Zoltan is compiled with 64-bit IDs. Since the resulting error message may be
      // hard to fathom, we warn the user ahead of time here.
#if defined(PHIST_HAVE_ISORROPIA)
      if (sizeof(phist_gidx)!=sizeof(int))
      {
        PHIST_SOUT(PHIST_WARNING,"WARNING: You are not using 32-bit integers as gids, Isorropia may run \n"
                                 "         into trouble and cause an obscure exception in Epetra.\n"
                                 "         To use 32-bit gids, pass -DPHIST_FORCE_INT_GIDX=ON to cmake.");
      }
#endif

#if 0

      //Now we're going to create a Epetra_Vector with weights to
      //be used as hypergraph edge weights in the repartitioning operation.
      //
      // We think of the rows as vertices of the hypergraph, and the columns
      // as hyperedges.  Our row matrix is square, so we can use the
      // the row map to indicate how the column weights should be
      // distributed.

      Teuchos::RCP<Epetra_Vector> hge_weights =
        Teuchos::rcp(new Epetra_Vector(rowmatrix->RowMatrixRowMap()));

      double* vals = hge_weights->Values();
      const Epetra_BlockMap& map = rowmatrix->RowMatrixRowMap();
      int num = map.NumMyElements();

      // we could put a weight according to the number of nonzeros in each row on the hyperedge
      for(int i=0; i<num; ++i) 
      {
        vals[i] = 1.0;
      }

      Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs;
      PHIST_TRY_CATCH(costs = Teuchos::rcp(new Isorropia::Epetra::CostDescriber()));

      PHIST_TRY_CATCH(costs->setHypergraphEdgeWeights(hge_weights),*iflag);

      //Now create the partitioner
      PHIST_TRY_CATCH(partitioner = Teuchos::rcp(new Isorropia::Epetra::Partitioner(rowmatrix, costs, paramlist)));
#else
      //Now create the partitioner without weights
      PHIST_TRY_CATCH(partitioner = Teuchos::rcp(new Isorropia::Epetra::Partitioner(rowmatrix, paramlist)),*iflag);
#endif    
      bal_map=Teuchos::null;
      PHIST_TRY_CATCH(bal_map=partitioner->createNewMap(),*iflag);

      // Results

    //  double bal0, bal1, cutn0, cutn1, cutl0, cutl1;

    //  int numProcs=rowmatrix->RowMatrixRowMap().Comm().NumProc();

      // Balance and cut quality before partitioning
    /*
      double goalWeight = 1.0 / (double)numProcs; 
      ispatest::compute_hypergraph_metrics(*rowmatrix, *costs, goalWeight,
                     bal0, cutn0, cutl0);
    */

      if (!redist) return;

      //Next create a Redistributor object and use it to create a repartitioned
      //copy of the matrix.

      Isorropia::Epetra::Redistributor rd(partitioner);

      //Use a try-catch block because Isorropia will throw an exception
      //if it encounters an error.

      PHIST_TRY_CATCH(bal_matrix = rd.redistribute(*rowmatrix),*iflag);

      return;
    }

  } // namespace epetra_internal
} // namespace phist
#endif // PHIST_HAVE_ISORROPIA
