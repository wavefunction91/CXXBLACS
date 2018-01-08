/*
 *  A simple C++ Wrapper for BLACS along with minimal extra functionality to 
 *  aid the the high-level development of distributed memory linear algebra.
 *  Copyright (C) 2016-2018 David Williams-Young

 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __INCLUDED_CXXBLACS_BLACSGRID_HPP__
#define __INCLUDED_CXXBLACS_BLACSGRID_HPP__

#include <cxxblacs/config.hpp>
#include <cxxblacs/blacs.hpp>
#include <cxxblacs/misc.hpp>
#include <cxxblacs/mpi.hpp>
#include <cxxblacs/lapack.hpp>

namespace CXXBLACS {


  struct LocalCoordinate {

    CB_INT locRowBlock; // Row Block of local buffer
    CB_INT locColBlock; // Col Block of local buffer

    CB_INT procRowOwn; // BLACS row of owning process
    CB_INT procColOwn; // BLACS col of owning process

    CB_INT locX; // Row coordinate within row block ( [0,MB) )
    CB_INT locY; // Col coordinate within col block ( [0,NB) )

  };

  /**
   * \brief A class to wrap the FORTRAN functionality of the BLACS (Basic 
   * Linear Algebra Communication Subroutines) and provide minimal extra 
   * functionality to aid in development.
   *
   * BlacsGrid is a class that is designed to wrap the FORTRAN routines of the 
   * BLACS (Basic Linear Algebra Communication Subroutines) that are used in 
   * conjunction with ScaLAPACK (Scalable Linear Albegra PACKage) and PBLAS 
   * (Parallel Basic Linear Algebra Subroutines). In addition to providing
   * simple C++ wrappers to the BLACS routines, BlacsGrid also provides extra 
   * (but minimal) functionality to aid the in the development of distributed 
   * memory linear algebra programs, such as:
   *   - Basic MPI execution design (such as ring execution)
   *   - Gathering a distributed matrix from the BLACS grid to root process
   *   - Distributing a matrix from a root process to the BLACS grid
   *   - Conversion of local and distrubuted memory element coordinates
   *
   * 
   */
  class BlacsGrid {

    // Global MPI information
    MPI_Comm comm_; ///< MPI Communicator
    CB_INT iProc_;  ///< Current MPI process (rank)
    CB_INT nProc_;  ///< Size of MPI universe (size)

    // BLACS grid information
    CB_INT IContxt_ = 0;  ///< BLACS context
    CB_INT nProcRow_;     ///< Number of rows on BLACS grid
    CB_INT nProcCol_;     ///< Number of columns on BLACS grid
    CB_INT iProcRow_;     ///< Current row coordinate of BLACS grid
    CB_INT iProcCol_;     ///< Current column coordinate of BLACS grid
    CB_INT mb_;           ///< Row block size (default = 2)
    CB_INT nb_;           ///< Column block size (default = 2)

  public:


    /**
     * \brief Constructor
     *
     * Initialize the BLACS grid. Attempts to create as close to a square
     * BLACS grid as possible. Obtains MPI and BLACS information about the
     * current proces, such as iProc, nProc, iProcRow, etc... The default
     * BLACS process grid is taken to be row-major
     *
     *   @param[in] c      MPI Communicator
     *   @param[in] MB     Block size for row distribution
     *   @param[in] NB     Block size for column distribution
     *   @param[in] ORDER  Process Grid ordering (row / column major)
     *
     */
    BlacsGrid(MPI_Comm c, CB_INT mb = 2, CB_INT nb = 2,
      std::string ORDER = "row-major") : 
      comm_(c), mb_(mb), nb_(nb) {

      // Check if MPI has been initialized
      int flag;
      MPI_Initialized(&flag);
      if( not flag ) {
        std::runtime_error err("MPI Environment Not Initialized!");
        throw err;
      }

      // Comm must be global for the time beging
      if( comm_ != MPI_COMM_WORLD ) {

        std::runtime_error err("BlacsGrid + non MPI_COMM_WORLD NYI");
        throw err;

      }

      BlacsPINFO(iProc_,nProc_); // Get global MPI info
      IContxt_ = BlacsGet(-1,0); // Setup BLACS context

      // Make as close to a square grid as possible
      nProcRow_ = int(std::sqrt(nProc_));
      nProcCol_ = nProc_ / nProcRow_;

      // Initialize BLACS grid
      BlacsGridInit(IContxt_,ORDER.c_str(),nProcRow_,nProcCol_);

      // Get grid information
      BlacsGridInfo(IContxt_,nProcRow_,nProcCol_,iProcRow_,iProcCol_);

    };




    // Getters 
    inline CB_INT iProc()    const { return iProc_;    }; ///< #iProc_
    inline CB_INT iProcRow() const { return iProcRow_; }; ///< #iProcRow_
    inline CB_INT iProcCol() const { return iProcCol_; }; ///< #iProcCol_
    inline CB_INT nProc()    const { return nProc_;    }; ///< #nProc_
    inline CB_INT nProcRow() const { return nProcRow_; }; ///< #nProcRow_
    inline CB_INT nProcCol() const { return nProcCol_; }; ///< #nProcCol_
    inline CB_INT iContxt()  const { return IContxt_;  }; ///< #IContxt_
    inline CB_INT NB()       const { return nb_;       }; ///< #nb_
    inline CB_INT MB()       const { return mb_;       }; ///< #mb_


    // Print functions
          
    // Print generic MPI / BLACS coordinate
    inline void printCoord(std::ostream &out) {

      std::stringstream ss;
      ss << "Greetings from process \n"; 
      ss << "  MPI:   " << iProc_ << " / " << nProc_ << "\n";
      ss << "  BLACS: ";
      ss << "(" << iProcRow_ << " , " << iProcCol_ << ") / ";
      ss << "(" << nProcRow_ << " , " << nProcCol_ << ")\n";
      out << ss.str();
      BlacsBarrier(IContxt_,"A");

    }


    // Print a local buffer in a user defined way
    template <typename Field, typename Op>
    inline void printLocalBuffer(const Op &op, const CB_INT M, 
      const CB_INT N, const Field *ALoc, const CB_INT LDLOCA) {

      CB_INT NLocR, NLocC;
      std::tie(NLocR,NLocC) = getLocalDims(M,N);
      
      RingExecute(comm_,[&](){ op(ALoc,NLocR,NLocC,LDLOCA)} );

    }

    // Print a local buffer in a standard way
    template <typename Field>
    inline void printLocalBuffer(std::ostream &out, const CB_INT M, 
      const CB_INT N, const Field *ALoc, const CB_INT LDLOCA) {

      CB_INT NLocR, NLocC;
      std::tie(NLocR,NLocC) = getLocalDims(M,N);
      
      RingExecute(comm_, [&]() {

        out << " Local Buffer (" << iProcRow_ << ", " << iProcCol_ << ")\n";
        out << std::scientific << std::setprecision(8);
        for(auto iLocR = 0; iLocR < NLocR; iLocR++) {
        for(auto iLocC = 0; iLocC < NLocC; iLocC++) {
          out << std::setw(20) << ALoc[iLocR + iLocC*LDLOCA] << ", ";
        }
        out << "\n";
        }
        out << "\n\n";

      });

    }


    // Misc helper functions
      
    inline INDX getLocalDims(const CB_INT M, const CB_INT N) const {

      return GetLocalDims(M,N,mb_,nb_,iProcRow_,iProcCol_,0,0,
          nProcRow_,nProcCol_);
             

    }

    inline INDX globalFromLocal(const CB_INT iX, const CB_INT iY) const {

      INDX indx;
      GlobalFromLocal(IContxt_,mb_,nb_,nProcRow_,nProcCol_,indx.first,
          indx.second,iProcRow_,iProcCol_,iX,iY);

      return indx;

    }

    inline LocalCoordinate localFromGlobal(const CB_INT I, const CB_INT J) 
      const {

      LocalCoordinate lc;
      LocalFromGlobal(IContxt_,mb_,nb_,nProcRow_,nProcCol_,I,J,
          lc.locRowBlock,lc.locColBlock,lc.procRowOwn,lc.procColOwn,
          lc.locX,lc.locY);

      return lc;

    }


    // Point-to-pont communication wrappers
      
    // Send
    template <typename Field>
    inline void Send(const CB_INT M, const CB_INT N, Field *A, const CB_INT LDA,
      const CB_INT DestRow, const CB_INT DestCol) const {

      GESD2D(IContxt_,M,N,A,LDA,DestRow,DestCol);

    }

    // Recv
    template <typename Field>
    inline void Recv(const CB_INT M, const CB_INT N, Field *A, const CB_INT LDA,
      const CB_INT OrigRow, const CB_INT OrigCol) const {

      GERV2D(IContxt_,M,N,A,LDA,OrigRow,OrigCol);

    }




    // Scatter / Gather

    template <typename Field>
    inline void Scatter( const CB_INT M, const CB_INT N, Field *A, 
      const CB_INT LDA, Field *ALoc, const CB_INT LDLOCA, 
      const CB_INT iSource, const CB_INT jSource ) {


      // If there's only one process, just copy the buffer
      if( nProcRow_ == 1 and nProcCol_ == 1 ) {

        LACOPY('A',M,N,A,LDA,ALoc,LDLOCA);
        return;

      }

      // Get local dimensions for current process
      CB_INT NLocR, NLocC;
      std::tie(NLocR,NLocC) = getLocalDims(M,N);


      // Source process block
      if( iProcRow_ == iSource and iProcCol_ == jSource ) {

        CB_INT I,J; // Global indicies

        for(CB_INT iPc = 0; iPc < nProcCol_; iPc++)
        for(CB_INT iPr = 0; iPr < nProcRow_; iPr++) {

          if( iPr == iSource and iPc == jSource ) continue;

          // Get local leading dimension
          CB_INT ldLocA;
          Recv(1,1,&ldLocA,1,iPr,iPc);

          // Get local variables for BLACS coordinate (iPr,iPc)
          CB_INT nLocR, nLocC;
          std::tie(nLocR,nLocC) = 
            GetLocalDims(M,N,mb_,nb_,iPr,iPc,0,0,nProcRow_,nProcCol_);
          

          // Check (non-robustly) if the local buffer is large enough
          // to fit the sub block of the matrix
          if( ldLocA > LDLOCA ) {

            std::stringstream ss;
            ss << "Local buffer does not seem to be large enough "
               << "to fit sub matrix in Scatter (" << iPr << ", " << iPc
               << " )";
            std::runtime_error err(ss.str());
              
            throw err;

          }

          // Copy local parts of the matrix to ALoc (temp buffer)
          for(CB_INT iLocC = 0; iLocC < nLocC; iLocC++)
          for(CB_INT iLocR = 0; iLocR < nLocR; iLocR++) {

            GlobalFromLocal(IContxt_,mb_,nb_,nProcRow_,nProcCol_,I,J,
              iPr,iPc,iLocR,iLocC);

            ALoc[iLocR + iLocC*ldLocA] = A[I + J*LDA];

          }

          // Send buffer to (iPr, iPc)
          Send(nLocR,nLocC,ALoc,ldLocA,iPr,iPc);

        } // Loop over other processes

        // Handle the local buffer
        for(CB_INT iLocC = 0; iLocC < NLocC; iLocC++)
        for(CB_INT iLocR = 0; iLocR < NLocR; iLocR++) {
          std::tie(I,J) = globalFromLocal(iLocR,iLocC);
          ALoc[iLocR + iLocC*LDLOCA] = A[I + J*LDA];
        }

      } else { // end Source process block

        // Non source process code
        Send(1,1,&LDLOCA,1,iSource,jSource); // Send LDLOCA to source
        Recv(NLocR,NLocC,ALoc,LDLOCA,iSource,jSource); // Recv buffer

      }


    } // Scatter




    template <typename Field>
    inline void Gather( const CB_INT M, const CB_INT N, Field *A, 
      const CB_INT LDA, Field *ALoc, const CB_INT LDLOCA, 
      const CB_INT iDest, const CB_INT jDest ) {


      // If there's only one process, just copy the buffer
      if( nProcRow_ == 1 and nProcCol_ == 1 ) {

        LACOPY('A',M,N,ALoc,LDLOCA,A,LDA);
        return;

      }

      // Get local dimensions for current process
      CB_INT NLocR, NLocC;
      std::tie(NLocR,NLocC) = getLocalDims(M,N);


      // Dest process block
      if( iProcRow_ == iDest and iProcCol_ == jDest ) {

        CB_INT I,J; // Global indicies


        // Handle the local buffer
        for(CB_INT iLocC = 0; iLocC < NLocC; iLocC++)
        for(CB_INT iLocR = 0; iLocR < NLocR; iLocR++) {
          std::tie(I,J) = globalFromLocal(iLocR,iLocC);
          A[I + J*LDA] = ALoc[iLocR + iLocC*LDLOCA];
        }




        // Loop over other processes
        for(CB_INT iPc = 0; iPc < nProcCol_; iPc++)
        for(CB_INT iPr = 0; iPr < nProcRow_; iPr++) {

          if( iPr == iDest and iPc == jDest ) continue;

          // Get local leading dimension
          CB_INT ldLocA;
          Recv(1,1,&ldLocA,1,iPr,iPc);

          // Get local variables for BLACS coordinate (iPr,iPc)
          CB_INT nLocR, nLocC;
          std::tie(nLocR,nLocC) = 
            GetLocalDims(M,N,mb_,nb_,iPr,iPc,0,0,nProcRow_,nProcCol_);
          

          // Check (non-robustly) if the local buffer is large enough
          // to fit the sub block of the matrix
          if( ldLocA > LDLOCA ) {

            std::stringstream ss;
            ss << "Local buffer does not seem to be large enough "
               << "to fit sub matrix in Gather (" << iPr << ", " << iPc
               << " )";
            std::runtime_error err(ss.str());
              
            throw err;

          }


          // Recv buffer from (iPr, iPc)
          Recv(nLocR,nLocC,ALoc,ldLocA,iPr,iPc);

          // Copy local parts of the matrix to ALoc (temp buffer)
          for(CB_INT iLocC = 0; iLocC < nLocC; iLocC++)
          for(CB_INT iLocR = 0; iLocR < nLocR; iLocR++) {

            GlobalFromLocal(IContxt_,mb_,nb_,nProcRow_,nProcCol_,I,J,
              iPr,iPc,iLocR,iLocC);

            A[I + J*LDA] = ALoc[iLocR + iLocC*ldLocA];

          }
        }

      } else { // end Dest process block

        // Non Dest process code
        Send(1,1,&LDLOCA,1,iDest,jDest); // Send LDLOCA to Dest
        Send(NLocR,NLocC,ALoc,LDLOCA,iDest,jDest); // Send buffer

      }


    } // Gather




  };

}; // CXXBLACS

#endif

