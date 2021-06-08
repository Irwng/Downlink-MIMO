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

/**********************************
 * basic constant model parameters
 **********************************/

constexpr double power = 1;
constexpr int M = 1;                              /* Number of time slots resources */
constexpr int Nt = 2;                             /* number of antennas at transmitter */
constexpr int Nr = 1;                             /* number of antennas at recevier */
constexpr int Mod = 2;				              /* BPSK modulation order */
constexpr double PI = 3.141592653589793;
    
constexpr int MinSNRdB = 0;
#ifdef DebugMode
    constexpr long NLoop = pow(10, 1);            /* number of simulation loops  */
    constexpr int MaxSNRdB = MinSNRdB;
#else
    constexpr long NLoop = pow(10, 6);            /* number of simulation loops  */
    constexpr long MaxSNRdB = 30;           
#endif

constexpr int Step = 3;

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* source codewords, 1*M */
typedef Matrix<int, 1, M> SourceMatrix;
extern SourceMatrix Source;
extern SourceMatrix Decode;

/* symbols after modulation, 1*M */
typedef Matrix<ComplexD, 1, M> ModuMatrix;
extern ModuMatrix Modu;

/* symbols after Alamouti-ST, Nt*M  */
typedef Matrix<ComplexD, Nt, M> STMatrix;
extern STMatrix SymAfterST;

/* channel parameters , Nr*Nt */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H;

/* signals after fadding channal, Nr*M */
typedef Matrix<ComplexD, Nr, M> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC;

/* signals after post processing, Nr*M */
/* After the post processing, the dimensionality of the signal matrix turns to Nt*M */
typedef Matrix<ComplexD, 1, M> SymAfterPPMatrix; 
extern SymAfterPPMatrix SymAfterPP;

/* signals after MAP, 1*M */
typedef Matrix<double, 1, M> SymAfterMAPMatrix;
extern SymAfterMAPMatrix SymAfterMAP;

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
 * description: Alamouti
 * date: 2021/3/16
 * input parameters: ModuMatrix& modu
 * output parameters: STMatrix& st
 ***************************************/
void Diversity(ModuMatrix& modu, STMatrix& st);


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
void Receiver(STMatrix& symAfterST,
              SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h,
              SymAfterPPMatrix& symAfterPP,
              SourceMatrix& source,  
              SourceMatrix& decode);

#endif //UDM_NOMA_DL_CLST_TEQUALR_V4_HEADER_H
