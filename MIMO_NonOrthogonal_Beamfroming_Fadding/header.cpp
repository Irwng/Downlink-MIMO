/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"
double N_Var;
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;                                     /* Bits Error Rate */
double power;
double PowerRate = 2;                               /* power allocation rate in beamforming matrix */
double beta; 
double sum = 0; 
fstream outfile;

SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
ModuMatrix MasterConstell[Mpoint];
CSIMatrix H[U];                                        
BFMatrix W;                                         /* beamforming matrix */
SymAfterFCMatrix SymAfterFC[U];                     /* receiving signals after postprocessing */
SymAfterPPMatrix SymAfterPP[U];                     /* receiving signals after preprocessing, Nt*1 */
SymAfterBFMatrix SymAfterBF;                        /* receiving signals after postprocessing */
DecodeMatrix Decode[U];  


void Initialize(ModuMatrix* masterConstell){

    outfile.open("Downlink MU-MIMO Beamforming.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    cout<<"Downlink MU-MIMO Beamforming"<<endl;
    outfile<<"Downlink MU-MIMO Beamforming"<<endl;

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;
    
#ifdef BPSK
    cout<<"BPSK"<<endl;   
    outfile<<"BPSK"<<endl; 
#endif
#ifdef QPSK 
    cout<<"QPSK"<<endl;   
    outfile<<"QPSK"<<endl;			            
#endif

#ifndef MAP
    cout<<"PowerRate: "<<PowerRate<<endl;
    outfile<<"PowerRate: "<<PowerRate<<endl;
#else
    cout<<"MAP"<<endl;
    outfile<<"MAP"<<endl;
#endif

    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl; 
   
    cout<<"SNRdB"<<setw(15)<<"BER"<<endl;
    outfile<<"SNRdB"<<setw(15)<<"BER"<<endl;

    /* initialize the Masterconstellation */
    /* maybe can simplified by bitset instead of this complex translation */
    for(int i = 0; i < Mpoint; ++i){                                     
        bitset<Nt> bit(i);
        for(int nt = 0; nt < Nt; ++nt){
            masterConstell[i](nt) = static_cast<double>(2 * bit[nt] - 1);
        }
    } 
}


void ChannelInitialize(int snrdB){

    /* the total power of transmitter is fixed */
    double snr = pow(10,(double)snrdB/10);
    power = 1.0;
    N_Var = power*Nt / snr;
    sum = 0.0;
    BER_TOTAL = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int i = 0; i < LenBit; i++){
        source(i)= rand()%2;
    }
    #ifdef DebugMode
        cout<<"source: "<<endl<<source<<endl;
    #endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    #ifdef BPSKMode
        /* BPSK: 0->-1, 1->1 */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD(2 * source(nt) - 1, 0);
        }
    #else
        /* QPSK: 00->a+aj, 01->-a+aj, 11-> -a-aj, 10->a-aj */
        /* QPSK: forward bits->1st antenna, next bits->2st antenna */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD((1 - 2 * source(BitperSymbol * nt + 1)) * sqrt(0.5), (1 - 2 * source(BitperSymbol * nt)) * sqrt(0.5));
        }
    #endif

    #ifdef DebugMode
        cout<<"modu"<<endl<<modu<<endl;
    #endif
}


void FadingChannel(CSIMatrix* h){

    /* channel parameters */
    for(int u = 0; u < U; ++u){
        for(int nr = 0; nr < Nr; ++nr){
            for(int nt = 0; nt < Nt; ++nt){
                double GG = sqrt(-2*log(randN() + 0.000000001));
                double B = randN();
                h[u](nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                        sqrt(0.5) * GG * sin(2*PI*B));
            }
        }
    }
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar/2.0)* GG * cos(2*PI*B), sqrt(nvar/2.0)* GG * sin(2*PI*B));
    return GuassN;
}


void Beamforming(ModuMatrix& modu, CSIMatrix* h, BFMatrix& w, 
                 SymAfterBFMatrix& symAfterBF, char* argv[]){

    double base = 0.0;
    
    /* initialize of W */ 
    for(int nt = 0; nt < Nt; ++nt){
        for(int nr = 0; nr < Nt; ++nr){
            w(nt, nr) = ComplexD(pow(PowerRate, nr), 0);
        }
    }

#ifdef DebugMode
    cout<<"w before normalization: "<<endl<<w<<endl;
#endif

    BFMatrix Denomiator_W;
    Denomiator_W = w * (w.conjugate().transpose());
    base = 0;
    for(int nt = 0; nt < Nt; ++nt)
        base += abs(Denomiator_W(nt,nt));
    beta = sqrt(power/base);
    w *= beta;

#ifdef DebugMode
    cout<<"w after normalization: "<<endl<<w<<endl;
#endif

    /* beamforming */
    symAfterBF = w * modu;

#ifdef DebugMode
    cout<<"symAfterBF"<<endl<<symAfterBF<<endl;
    sum = 0;
    for(int nt = 0; nt < Nt; ++nt){
        sum += pow(abs(symAfterBF(nt)), 2);
    }
    cout<<"weight of symAfterBF: "<<sum<<endl;
#endif
}


void Receiver(SymAfterBFMatrix& symAfterBF, SymAfterFCMatrix* symAfterFC, CSIMatrix* h, 
              SymAfterPPMatrix* symAfterPP, BFMatrix& w, SourceMatrix& source,
              DecodeMatrix* decode, ModuMatrix* masterConstell){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; ++nr){
        tmp(nr) = AWGN(N_Var);
    }
    /* post-processing */
    ComplexD gamma;
    for(int u = 0; u < U; ++u){
        symAfterFC[u] = h[u] * symAfterBF + tmp;// 

        /* 1st step: processing the signal */
        ComplexD CSISum = p0;
        for(int nt = 0; nt < Nt; ++nt)
            CSISum += h[u](0, nt);
        symAfterFC[u] /= (CSISum*beta);
    
        /* 2nd step: detect the signal with bigger power */
        symAfterPP[u](1) = symAfterFC[u](0) / PowerRate;
    #ifdef BPSKMode
        decode[u](1) = static_cast<int>(symAfterPP[u](1).real() > 0);
        BER_TOTAL += (decode[u](1)!=source(1));

    #else
        decode[u](BitperSymbol) = static_cast<int>(symAfterPP[u](1).imag() < 0);
        BER_TOTAL += (decode[u](BitperSymbol)!=source(BitperSymbol));
        decode[u](BitperSymbol + 1) = static_cast<int>(symAfterPP[u](1).real() < 0);
        BER_TOTAL += (decode[u](BitperSymbol + 1)!=source(BitperSymbol + 1));
    #endif

        /* 3rd step: do SIC to detect the rest signal */
        symAfterPP[u](0) = symAfterFC[u](0) - PowerRate * ComplexD(2 * decode[u](1) - 1, 0);
    #ifdef BPSKMode
        decode[u](0) = static_cast<int>(symAfterPP[u](0).real() > 0);
        BER_TOTAL += (decode[u](0)!=source(0));
    #else
        decode[u](0) = static_cast<int>(symAfterPP[u](0).imag() < 0);
        BER_TOTAL += (decode[u](0)!=source(0));
        decode[u](1) = static_cast<int>(symAfterPP[u](0).real() < 0);
        BER_TOTAL += (decode[u](1)!=source(1));
    #endif

    #ifdef DebugMode
        if(u == 0){
            cout<<"data for user0"<<endl;
            cout<<"origin SymAfterFC[0]"<<endl<<h[0] * symAfterBF<<endl;
            cout<<"h[0]"<<endl<<h[0]<<endl;
            cout<<"CSISum"<<endl<<CSISum<<endl;
            cout<<"SymAfterFC[0] after post processing"<<endl<<symAfterFC[0]<<endl;
            cout<<"SymAfterPP[0]"<<endl<<symAfterPP[0]<<endl;
            cout<<"Decode[0]"<<endl<<decode[0]<<endl;
        }
    #endif
    }
}


// #else
    //     SymAfterFCMatrix MasterConstellFixed[Mpoint];
    // 	double d[Mpoint];                        /* Euclidean dist of each star point*/  
    //     double MinDist = INT16_MAX;                    /* the minimum Euclidean dist */
    //     int MinDistPoint = 0;                    /* the minimum Euclidean dist point */

    //     for(int i = 0; i < Mpoint; ++i){                                     
    //         MasterConstellFixed[i] = h[u] * masterConstell[i];
    //         d[i] = 0;
    //         for(int nr = 0; nr < Nr; ++nr){
    //             d[i] += pow(abs(symAfterFC[u](nr) - MasterConstellFixed[i](nr)), 2);
    //         }
    //         if(d[i] < MinDist){
    //             MinDist = d[i];
    //             MinDistPoint = i;
    //         }
    //     }

    //     /* completely correct */
    //     bitset<Nt> bits(MinDistPoint);

    //     for(int nt = 0; nt < Nt; ++nt){
    //         decode[u](nt) = bits[nt];
    //         BER_TOTAL += (decode[u](nt)!=source(nt));
    //     }
    // #ifdef DebugMode
    //     if(u == 0){
    //         cout<<"origin SymAfterFC[0]"<<endl<<h[0] * symAfterBF<<endl;
    //         cout<<"MasterConstellFixed"<<endl;
    //         for(int i = 0; i < Mpoint; ++i) cout<<MasterConstellFixed[i]<<endl;
    //         cout<<"MinDistPoint"<<endl<<MinDistPoint<<endl;
    //         cout<<"bits"<<endl<<bits<<endl;
    //         cout<<"decode[u]"<<endl<<decode[0]<<endl;
    //     }
    // #endif

    // #endif