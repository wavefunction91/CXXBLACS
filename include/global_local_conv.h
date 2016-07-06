// Member functions
/**
 * \brief Calculate local 2D coordinates from global (distributed) 2D
 *        coordinates using the current processes BLACS coordinates
 *
 * @param[in]  I Distributed row coordinate
 * @param[in]  J Distributed column coordinate
 * @param[out] iX Local row coordinate
 * @param[out] iY Local column coordinate
 */
inline void localFromGlobal(const int I, const int J, int &iX, int &iY) {

  int L,M;
  int tmpPc(this->iProcCol_), tmpPr(this->iProcRow_);
  BlacsGrid::LocalFromGlobal(this->IContxt_,this->mb_,this->nb_,
    this->nProcRow_,this->nProcCol_,I,J,L,M,tmpPr,tmpPc,iX,iY);

  iX = L*this->mb_ + iX;
  iY = M*this->nb_ + iY;
}

/**
 * \brief Calculate global (distributed) 2D coordinates from local 2D
 *        coordinates using the current process' BLACS coordinates
 *
 * @param[out]  I Distributed row coordinate
 * @param[out]  J Distributed column coordinate
 * @param[in]   iX Local row coordinate
 * @param[in]   iY Local column coordinate
 */
inline void globalFromLocal(int &I, int &J, const int iX, const int iY) {

  int L,M;
  BlacsGrid::GlobalFromLocal(this->IContxt_,this->mb_,this->nb_,
    this->nProcRow_,this->nProcCol_,I,J,this->iProcRow_,this->iProcCol_,iX,iY);

}

/**
 * \brief Obtain dimension of local buffer of a distributed matrix using current
 * process' BLACS coordinates
 *
 * @param [in]  N Number of rows of distributed matrix
 * @param [in]  M Number of columns of distributed matrix
 * @param [out] NLocR Number of rows of local buffer
 * @param [out] NLocC Number of columns of local buffer
 */
inline void getLocalDims(const int N, const int M, int &NLocR, 
  int &NLocC) {
  BlacsGrid::GetLocalDims(N,M,this->mb_,this->nb_,this->nProcRow_,
    this->nProcCol_,this->iProcRow_,this->iProcCol_,NLocR,NLocC);
}

// Helper functions
/**
 * \brief Calculate local 2D coordinates from global (distributed) 2D
 *        coordinates using general BLACS grid and coordinates
 *
 * @param[in]  ICONTXT BLACS context
 * @param[in]  MB Distributed row block size
 * @param[in]  NB Distributed column block size
 * @param[in]  NPROW Number of rows in process grid
 * @param[in]  NPCOL Number of columnss in process grid
 * @param[in]  I Distributed row coordinate
 * @param[in]  J Distributed column coordinate
 * @param[out] L Local row block of local buffer
 * @param[out] M Local column block of local buffer
 * @param[out] Pr Process row the contains distributed element (I,J)
 * @param[out] Pc Process column the contains distributed element (I,J)
 * @param[out] iX Local row coordinate
 * @param[out] iY Local column coordinate
 */
static inline void LocalFromGlobal(const int ICONTXT, const int MB, 
  const int NB, const int NPROW, const int NPCOL, const int I, const int J, 
  int &L, int &M, int &Pr, int &Pc, int &iX, int &iY) {
  
  L  = I / (NPROW * MB);
  M  = J / (NPCOL * NB);
  Pr = (I / MB) % NPROW;
  Pc = (J / NB) % NPCOL;
  iX = I % MB;
  iY = J % NB;
}

/**
 * \brief Calculate global (distributed) 2D coordinates from local
 *        coordinates using general BLACS grid and coordinates
 *
 * @param[in]  ICONTXT BLACS context
 * @param[in]  MB Distributed row block size
 * @param[in]  NB Distributed column block size
 * @param[in]  NPROW Number of rows in process grid
 * @param[in]  NPCOL Number of columnss in process grid
 * @param[out] I Distributed row coordinate
 * @param[out] J Distributed column coordinate
 * @param[in]  Pr Process row the contains distributed element (I,J)
 * @param[in]  Pc Process column the contains distributed element (I,J)
 * @param[in]  iX Local row coordinate
 * @param[in]  iY Local column coordinate
 */
static inline void GlobalFromLocal(const int ICONTXT, const int MB,
  const int NB, const int NPROW, const int NPCOL, int &I, int &J,
  const int &Pr, const int &Pc, const int &iX, const int &iY) {

  int L = iX / MB;
  int M = iY / NB;

  I = (L * (NPROW - 1) + Pr) * MB + iX;
  J = (M * (NPCOL - 1) + Pc) * NB + iY;
}

