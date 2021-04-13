/***********************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the declaration of constants、variables、functions，
             and the definition of new type
************************************************************************/
#ifndef UDM_DL_CLST_TEQUALR_MPA_VITERBI_V3_HEADER_H
#define UDM_DL_CLST_TEQUALR_MPA_VITERBI_V3_HEADER_H

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

#define randN() (rand()/(double)RAND_MAX)         /* Random value in [0,30000] */
typedef complex<double> ComplexD;

/**********************************************************************************
 * Control Center: must choose one and only one of the choices in different blocks
 **********************************************************************************/

/* running mode*/
// #define DebugMode
#define MonteCarlo

/* channel */
// #define FaddingChannel
#define AWGNMode

/**********************************
 * basic constant model parameters
 **********************************/

constexpr int Nt = 1;                             /* number of antennas at transmitter */
constexpr int Nr = 1;                             /* number of antennas at recevier */
constexpr int Mod = 2;				              /* BPSK modulation order */
constexpr int Mpoint = pow(Mod, Nt);
constexpr double PI = 3.141592653589793;
constexpr double power = 4;
    
constexpr int minEbN0dB = 0;
#ifdef DebugMode
    constexpr long NLoop = pow(10, 0);            /* number of simulation loops  */
    constexpr int maxEbN0dB = minEbN0dB;
#else
    constexpr long NLoop = pow(10, 6);            /* number of simulation loops  */
    constexpr long maxEbN0dB = 30;           
#endif
constexpr int step = 3;

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* source codewords, 1*LenBit*/
typedef Matrix<int, Nt, 1> SourceMatrix;
extern SourceMatrix Source;
extern SourceMatrix Decode;

/* symbols after modulation, Nt*1  */
typedef Matrix<ComplexD, Nt, 1> ModuMatrix;
extern ModuMatrix Modu;
extern ModuMatrix Constell[Mpoint];

/* channel parameters , Nr*Nt */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H;

/* signals after fadding channal, Nr*1 */
typedef Matrix<ComplexD, Nr, 1> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC;
extern SymAfterFCMatrix ConstellFixed[Mpoint];

/*************************************
 * basic global variable declaration
 *************************************/

extern double N_Var;                      /* variance of Noise*/
extern double BER_TOTAL;                  /* total number of error bits*/
extern double BER;                        /* total number of error bits*/
extern fstream outfile;

/***********************************************
 * basic functions declaration
 ***********************************************/

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void NormalIO();
void InitMapMatrix();
/**************************************
 * description: add AAWGN noise
 * date: 2020/8/20
 * input parameters: variance of noise
 * output parameters: AWGN noise in complex
 ***************************************/
ComplexD AWGN(double nvar);


/**************************************
 * description: overload the operator - between ComplexD & int
 * date: 2020/11/20
 * input parameters: ComplexD comp, int b
 * output parameters: ComplexD tmp
 ***************************************/
ComplexD operator-(ComplexD comp, int b);


/**************************************
 * description: calculate the Eb/N0 & initialize the error counter
 * date: 2020/8/18
 * input parameters: ebN0dB
 * output parameters: N_Var
 ***************************************/
void ChannelInitialize(int ebN0dB);


/**************************************
 * description: generater the source bits
 * date: 2020/9/24
 * input parameters: number of the users(U)
 * output parameters: the source bits(Source[U])
 ***************************************/
void BitSource(SourceMatrix& source);


/**************************************
 * description: Modulation
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
 * description: receiving signals
 * date: 2020/8/22
 * input parameters: SymAfterBFMatrix* symafterbf
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Receiver(ModuMatrix& modu, 
              CSIMatrix& h,
              ModuMatrix* constell,
              SourceMatrix& source,  
              SourceMatrix& decode);

#endif //UDM_NOMA_DL_CLST_TEQUALR_V4_HEADER_H
