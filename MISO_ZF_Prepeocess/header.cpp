/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"

double power;                                       
double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;                                     /* Bits Error Rate */
fstream outfile;

SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
CSIMatrix H;                                        
ProcessMatrix V;                                        /* preprocess matrix */
SymAfterFCMatrix SymAfterFC;                        /* receiving signals */ 
ModuMatrix SymAfterPP;                        /* receiving signals after postprocessing */
DecodeMatrix Decode;  

void Initialize(){

    outfile.open("Downlink MIMO Preprocess.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    cout<<"Downlink MIMO Preprocess"<<endl;
    outfile<<"Downlink MIMO Preprocess"<<endl;

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;

    cout<<"BPSK"<<endl;   
    outfile<<"BPSK"<<endl; 

    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"NLoop: "<<NLoop<<endl; 
   
    cout<<"SNRdB"<<setw(15)<<"BER"<<endl;
    outfile<<"SNRdB"<<setw(15)<<"BER"<<endl;
    
}


void ChannelInitialize(int snrdB){

    /* the total power of transmitter is fixed */
    double snr = pow(10,(double)snrdB/10);
    power = 1.0;
    N_Var = Nt * power / snr;
    BER_TOTAL = 0;

#ifdef DebugMode
    cout<<"N_Var: "<<N_Var<<endl;
#endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int nt = 0; nt < Nt; ++nt){
        source(nt)= rand()%2;
    }
#ifdef DebugMode
    cout<<"source: "<<endl<<source<<endl;
#endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int nt = 0; nt < Nt; ++nt){
        modu(nt) = ComplexD(2 * source(nt) - 1, 0);
    }
#ifdef DebugMode
    cout<<"modulation: "<<endl<<modu<<endl;
#endif  

}


void FadingChannel(CSIMatrix& h){

    /* channel parameters */
    for(int nr = 0; nr < Nr; ++nr){
        for(int nt = 0; nt < Nt; ++nt){
            double GG = sqrt(-2*log(randN() + 0.000000001));
            double B = randN();
            h(nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                sqrt(0.5) * GG * sin(2*PI*B));
        }
    }

#ifdef DebugMode
    cout<<"h"<<endl<<h<<endl;
#endif 
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar/2)* GG * cos(2*PI*B), sqrt(nvar/2)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(ModuMatrix& modu, SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h, ProcessMatrix& v, SourceMatrix& source, 
              DecodeMatrix& decode, ModuMatrix& symAfterPP){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; ++nr){
        tmp(nr) = AWGN(N_Var);
    }
    symAfterFC = h * modu + tmp;
    v(0) = conj(h(0))/(h(0)*conj(h(0))+N_Var);
    v(1) = conj(h(1))/(h(1)*conj(h(1))+N_Var);
    symAfterPP = v * symAfterFC;

    /* check the performance of the post processing */
    for(int nt = 0; nt < Nt; ++nt){
        decode(nt) = static_cast<int>(symAfterPP(nt).real() > 0);
        BER_TOTAL += static_cast<double>(decode(nt) != source(nt));
    }
}
