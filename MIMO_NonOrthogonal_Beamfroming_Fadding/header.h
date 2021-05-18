/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.h
***********************************************************/
#ifndef MIMO_BEAMFORMING_FADDING_HEADER_H
#define MIMO_BEAMFORMING_FADDING_HEADER_H

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

#define randN() (rand()/(double)RAND_MAX)         /* Random value in [0, 1] */
typedef complex<double> ComplexD;

/* running mode*/
// #define DebugMode
#define MonteCarlo

/* Modulation mode*/
#define BPSKMode
// #define QPSKMode

/* MAP mode*/
// #define MAP

/* basic model parameters */ 

constexpr int U = 6;                              /* Number of users */
constexpr int Nr = 1;                             /* number of antennas at recevier */
constexpr int Nt = 2;                             /* number of antennas at transmitter */    
constexpr int LenBit = Nt;                        /* number of bits of all users */
constexpr int Mpoint = pow(2, Nt);                        /* Master-constellation Points */

#ifdef BPSKMode
    constexpr int Mod = 2;				          /* BPSK modulation order */
    constexpr int BitperSymbol = 1;			
#else
    constexpr int Mod = 4;				          /* QPSK modulation order */
    constexpr int BitperSymbol = 2;				            
#endif

#ifdef DebugMode
    constexpr int MinSNRdB = 30;
    constexpr long NLoop = pow(10, 0);            /* number of simulation loops  */
    constexpr int MaxSNRdB = MinSNRdB;
#else
    constexpr int MinSNRdB = 0;
    constexpr long NLoop = pow(10, 6);            /* number of simulation loops  */
    constexpr int MaxSNRdB = 30;           
#endif

constexpr int Step = 5;           
constexpr double PI = 3.141592653589793;
constexpr ComplexD p0(0, 0);

extern double N_Var;
extern double BER_TOTAL;                          /* total number of error bits*/
extern double BER;                                /* total number of error bits*/
extern double power;
extern double PowerRate;                          /* power allocation rate in beamforming matrix */
extern double beta;                               
extern fstream outfile;

/* source codewords, Nt*1 */
typedef Matrix<int, Nt * BitperSymbol, 1> SourceMatrix;
extern SourceMatrix Source;

/* symbols after modulation, Nt*1 */
typedef Matrix<ComplexD, Nt, 1> ModuMatrix;
extern ModuMatrix Modu;
extern ModuMatrix MasterConstell[Mpoint];

/* channel parameters , Nr*Nt */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H[U];

/* beamforming matrix, Nt*Nt */
typedef Matrix<ComplexD, Nt, Nt> BFMatrix;
extern BFMatrix W;

/* signals after Beamforming, Nt*1 */
typedef ModuMatrix SymAfterBFMatrix;
extern SymAfterBFMatrix SymAfterBF;

/* signals after fadding channal, Nr*1 */
typedef Matrix<ComplexD, Nr, 1> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC[U];

/* signals after preprocessing, Nt*1 */
typedef ModuMatrix SymAfterPPMatrix; 
extern SymAfterPPMatrix SymAfterPP[U];

/* final decoded results */
typedef SourceMatrix DecodeMatrix;
extern DecodeMatrix Decode[U];


/**************************************
 * description: initialize the io & parameters
 * date: 201/3/23
 * input parameters: char* Argv[]
 * output parameters: void
 ***************************************/
void Initialize(ModuMatrix* masterConstell);


/**************************************
 * description: calculate the Eb/N0 & initialize the error counter
 * date: 2020/8/18
 * input parameters: SNR
 * output parameters: Eb/N0
 ***************************************/
void ChannelInitialize(int snrdB);


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
void FadingChannel(CSIMatrix* h);


/**************************************
 * description: add AAWGN noise
 * date: 2020/8/20
 * input parameters: variance of noise
 * output parameters: AWGN noise in complex
 ***************************************/
ComplexD AWGN(double nvar);


/**************************************
 * description: receiving signals
 * date: 2021/3/17
 * input parameters: SymAfterBFMatrix* symafterbf
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Beamforming(ModuMatrix& modu,
                 CSIMatrix* h,
                 BFMatrix& w,
                 SymAfterBFMatrix& symAfterBF,
                 char* argv[]);


/**************************************
 * description: receiving signals
 * date: 2020/8/22
 * input parameters: SymAfterBFMatrix* symafterbf
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Receiver(SymAfterBFMatrix& symAfterBF,
              SymAfterFCMatrix* symAfterFC, 
              CSIMatrix* h,
              SymAfterPPMatrix* symAfterPP, 
              BFMatrix& w,
              SourceMatrix& source,
              DecodeMatrix* decode,
              ModuMatrix* masterConstell);

#endif //MIMO_BEAMFORMING_FADDING_HEADER_H
