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
#include <cxxblacs/scalapack.hpp>

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
    CB_INT bHandle_ = 0;
    CB_INT nProcRow_;     ///< Number of rows on BLACS grid
    CB_INT nProcCol_;     ///< Number of columns on BLACS grid
    CB_INT iProcRow_;     ///< Current row coordinate of BLACS grid
    CB_INT iProcCol_;     ///< Current column coordinate of BLACS grid
    CB_INT mb_;           ///< Row block size (default = 2)
    CB_INT nb_;           ///< Column block size (default = 2)

    CB_INT iSrc_;         ///< Source Row
    CB_INT jSrc_;         ///< Source Col

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
    BlacsGrid(MPI_Comm c, CB_INT mb, CB_INT nb,
      CB_INT npr = 0, CB_INT npc = 0, 
      std::string ORDER = "row-major", CB_INT iSrc = 0,
      CB_INT jSrc = 0) : 
      comm_(c), mb_(mb), nb_(nb), iSrc_(iSrc), jSrc_(jSrc) {

      // Check if MPI has been initialized
      int flag;
      MPI_Initialized(&flag);
      if( not flag ) {
        std::runtime_error err("MPI Environment Not Initialized!");
        throw err;
      }

      if( comm_ == MPI_COMM_NULL ) return;

      // Get the MPI info and system context
      int IPROC, NPROC; // for 64-bit ints
      MPI_Comm_rank(comm_,&IPROC);
      MPI_Comm_size(comm_,&NPROC);
      iProc_ = IPROC;
      nProc_ = NPROC;
      bHandle_ = Csys2blacs_handle( comm_ );
      IContxt_ = bHandle_;

      // Make as close to a square grid as possible
      if( npr && npc ) {        

        if( npr * npc != nProc_ ) {
          std::runtime_error err("NPROCROW * NPROCCOL != NPROC");
          throw err;
        }

        nProcRow_ = npr;
        nProcCol_ = npc;

      } else if( not ORDER.compare("linear") ) {

        nProcRow_ = 1;
        nProcCol_ = nProc_;
        ORDER = "row-major";

      } else {

        nProcRow_ = int(std::sqrt(nProc_));
        nProcCol_ = nProc_ / nProcRow_;

        while( nProcRow_ * nProcCol_ != nProc_ ) {                                  
          
          nProcRow_--;
          nProcCol_ = nProc_ / nProcRow_;                                           
        
        }

      }

      // Initialize BLACS grid
      BlacsGridInit(IContxt_,ORDER.c_str(),nProcRow_,nProcCol_);

      // Get grid information
      BlacsGridInfo(IContxt_,nProcRow_,nProcCol_,iProcRow_,iProcCol_);

    };


    ~BlacsGrid() {  
      BlacsGridExit(IContxt_); 
      Cfree_blacs_system_handle( bHandle_ ); 
    }




    // Getters 
    inline CB_INT iProc()    const noexcept { return iProc_;    }; ///< #iProc_
    inline CB_INT iProcRow() const noexcept { return iProcRow_; }; ///< #iProcRow_
    inline CB_INT iProcCol() const noexcept { return iProcCol_; }; ///< #iProcCol_
    inline CB_INT nProc()    const noexcept { return nProc_;    }; ///< #nProc_
    inline CB_INT nProcRow() const noexcept { return nProcRow_; }; ///< #nProcRow_
    inline CB_INT nProcCol() const noexcept { return nProcCol_; }; ///< #nProcCol_
    inline CB_INT iContxt()  const noexcept { return IContxt_;  }; ///< #IContxt_
    inline CB_INT NB()       const noexcept { return nb_;       }; ///< #nb_
    inline CB_INT MB()       const noexcept { return mb_;       }; ///< #mb_

    inline bool i_participate() const noexcept { return comm_ != MPI_COMM_NULL; };

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
      
      RingExecute(comm_,[&](){ op(ALoc,NLocR,NLocC,LDLOCA); } );

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

      return GetLocalDims(M,N,mb_,nb_,iProcRow_,iProcCol_,iSrc_,jSrc_,
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

      auto a = const_cast<typename std::remove_cv<Field>::type *>(A);
      GESD2D(IContxt_,M,N,a,LDA,DestRow,DestCol);

    }

    // Recv
    template <typename Field>
    inline void Recv(const CB_INT M, const CB_INT N, Field *A, const CB_INT LDA,
      const CB_INT OrigRow, const CB_INT OrigCol) const {

      auto a = const_cast<typename std::remove_cv<Field>::type *>(A);
      GERV2D(IContxt_,M,N,a,LDA,OrigRow,OrigCol);

    }


    // Broadcast
    template <typename Field>
    inline void Broadcast(const char SCOPE[], const char TOP[], const CB_INT M,
      const CB_INT N, Field *A, const CB_INT LDA, const CB_INT ISRC, 
      const CB_INT JSRC) {

      if( iProcRow_ == ISRC and iProcCol_ == JSRC )
        GEBS2D(IContxt_,SCOPE,TOP,M,N,A,LDA);
      else
        GEBR2D(IContxt_,SCOPE,TOP,M,N,A,LDA,ISRC,JSRC);

    };


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

      // Create a new BLACS grid which conains matrix
      // on source process
      BlacsGrid new_grid( comm_, M, N, 0,0, "row-major", 
                          iSource, jSource );

      
      // Get Descriptors of A on both grids
      auto DescA = 
        new_grid.descInit( M, N, iSource, jSource, LDA );
      auto DescALoc = 
        this->descInit( M, N, iSrc_, jSrc_, LDLOCA );


      // Redistribute
      PGEMR2D( M, N, 
               A,    1, 1, DescA, 
               ALoc, 1, 1, DescALoc,
               this->iContxt() );
               
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

      // Create a new BLACS grid which conains matrix
      // on dest process
      BlacsGrid new_grid( comm_, M, N, 0,0, "row-major", 
                          iDest, jDest );


      // Get Descriptors of A on both grids
      auto DescA = 
        new_grid.descInit( M, N, iDest, jDest, LDA );
      auto DescALoc = 
        this->descInit( M, N, iSrc_, jSrc_, LDLOCA );

      // Redistribute
      PGEMR2D( M, N, 
               ALoc, 1, 1, DescALoc,
               A,    1, 1, DescA, 
               this->iContxt() );


    } // Gather





    // ScaLAPACK Helper functions
      
    inline ScaLAPACK_Desc_t descInit(const CB_INT M, const CB_INT N,
      const CB_INT ISRC, const CB_INT JSRC, const CB_INT LDD) {

      return DescInit(M,N,mb_,nb_,ISRC,JSRC,IContxt_,LDD);

    };

  };




}; // CXXBLACS

#endif

