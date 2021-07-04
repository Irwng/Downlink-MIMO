/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.h
***********************************************************/
#ifndef UDM_NOMA_DL_CLST_TEQUALR_V4_HEADER_H
#define UDM_NOMA_DL_CLST_TEQUALR_V4_HEADER_H


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <bitset>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define randN() (rand()/(double)RAND_MAX)     /* Random value in [0, 1] */

/* running mode*/
// #define DebugMode
#define MonteCarlo

typedef complex<double> ComplexD;

/* basic model parameters */ 
constexpr int Nr = 1;                           /* number of antennas at recevier */
constexpr int Nt = 2;                          /* number of antennas at transmitter */    
constexpr int Mod = 2;				        /* BPSK modulation order */
constexpr int BitperSymbol = 1;			
constexpr int Mpoint = pow(Mod, Nt);

constexpr double PI = 3.141592653589793;
constexpr int MinSNRdB = 0;
#ifdef DebugMode
    constexpr long NLoop = pow(10, 1);          /* number of simulation loops  */
    constexpr int MaxSNRdB = MinSNRdB;
#else
    constexpr long NLoop = pow(10, 5);          /* number of simulation loops  */
    constexpr int MaxSNRdB = 30;           
#endif

constexpr int Step = 3;           

extern double power;
extern double N_Var;                            /* variance of Noise*/
extern double BER_TOTAL;                        /* total number of error bits*/
extern double BER;                              /* total number of error bits*/
extern fstream outfile;

/* source codewords, Nt*1 */
typedef Matrix<int, Nt * BitperSymbol, 1> SourceMatrix;
extern SourceMatrix Source;

/* symbols after modulation, Nt*1 */
typedef Matrix<ComplexD, Nt, 1> ModuMatrix;
extern ModuMatrix Modu;
extern ModuMatrix SymAfterPP;

/* channel parameters , Nr*Nt */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H;
typedef Matrix<ComplexD, Nt, Nr> ProcessMatrix;
extern ProcessMatrix V;

/* signals after fadding channal, Nr*1 */
typedef Matrix<ComplexD, Nr, 1> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC;

/* final decoded results */
typedef SourceMatrix DecodeMatrix;
extern DecodeMatrix Decode;


/**************************************
 * description: initialize the io & parameters
 * date: 2021/3/23
 ***************************************/
void Initialize();


/**************************************
 * description: calculate the Eb/N0 & initialize the error counter
 * date: 2020/8/18
 * input parameters: SNR
 * output parameters: Eb/N0
 ***************************************/
void ChannelInitialize(int SNR);


/**************************************
 * description: generater the source bits
 * date: 2020/9/24
 * input parameters: number of the users(U)
 * output parameters: the source bits(Source[U])
 ***************************************/
void BitSource(SourceMatrix& source);


/**************************************
 * description: initialize the codebooks & constellations
 * date: 2020/10/18
 * input parameters: CodeMatrix& code
 * output parameters: ModuMatrix& modu
 ***************************************/
void Modulation(SourceMatrix& source, ModuMatrix& modu);


/**************************************
 * description: flat-fading channel
 * date: 2020/8/16
 * input parameters: transmitting signals
 * output parameters: receiving signals pointer
 ***************************************/
void FadingChannel(CSIMatrix& h);


/**************************************
 * description: add AAWGN noise
 * date: 2020/8/20
 * input parameters: variance of noise
 * output parameters: AWGN noise in complex
 ***************************************/
ComplexD AWGN(double nvar);


/**************************************
 * description: receiving signals
 * date: 2020/8/22
 * input parameters: SymAfterBFMatrix* symafterbf
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Receiver(ModuMatrix& modu,
              SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h,
              ProcessMatrix& v,
              SourceMatrix& source,
              DecodeMatrix& decode,
              ModuMatrix& symAfterPP);

#endif //UDM_NOMA_DL_CLST_TEQUALR_V4_HEADER_H
