/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var = 0;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;

SourceMatrix Source;                         /* source codewords */
ModuMatrix Modu;                             /* symbols after modulation */
STMatrix SymAfterST;
CSIMatrix H;                                 /* channel parameters , U*Nt*Nr */
SymAfterFCMatrix SymAfterFC;
SymAfterPPMatrix SymAfterPP;                 /* receiving signals after postprocessing, Nj*Nr*M */
SymAfterMAPMatrix SymAfterMAP;
SourceMatrix Decode;  


void NormalIO(){
    
    cout<<"MIMO-Diversity"<<endl;
    outfile.open("MIMO-Diversity.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"MIMO-Diversity"<<endl;
    
    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"T: "<<M<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"T: "<<M<<setw(10)<<"NLoop: "<<NLoop<<endl; 
    

    cout<<"FadingChannel"<<endl;
    outfile<<"FadingChannel"<<endl;

    cout<<"snrdB"<<setw(15)<<"BER"<<endl;
    outfile<<"snrdB"<<setw(15)<<"BER"<<endl;
}


/* key point */
void ChannelInitialize(int snrdB){

    double snr = pow(10, (double)snrdB/10);
    N_Var = Nt * power / snr; //
    BER_TOTAL = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int m = 0; m < M; ++m){
        source(m)= rand()%2;
    }

    // #ifdef DebugMode
    //     cout<<".................source................."<<endl;
    //     cout<<source<<endl;
    // #endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int m = 0; m < M; ++m){
        modu(m) = ComplexD(2*source(m)-1, 0);
    }
    #ifdef DebugMode
        cout<<".................bpsk signal................."<<endl;
        cout<<modu<<endl;
    #endif
}


void Diversity(ModuMatrix& modu, STMatrix& st){

    for(int nt = 0; nt < Nt; ++nt){
        for(int m = 0; m < M; ++m){
            st(nt, m) = modu(m);
        }
    }

    #ifdef DebugMode
        cout<<".................st................."<<endl;
        cout<<st<<endl;
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
        cout<<"................h.................. "<<endl;
        cout<<h<<endl; 
    #endif
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar/2)* GG * cos(2*PI*B), sqrt(nvar/2)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(STMatrix& symAfterST, SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h, SymAfterPPMatrix& symAfterPP, 
              SourceMatrix& source,  SourceMatrix& decode){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; nr++){
        for(int m = 0; m < M; m++){
            tmp(nr, m) = AWGN(N_Var);
        }
    }
    symAfterFC = h * symAfterST + tmp;//
    ComplexD h_sum[Nr]; 
    double h_conj_sum = 0;
    for(int nr = 0; nr < Nr; ++nr){
        h_sum[nr] = ComplexD(0,0);
        for(int nt = 0; nt < Nt; ++nt){
            h_sum[nr] += h(nr, nt);
        }
        h_conj_sum += pow(abs(h_sum[nr]),2);
    }

    /* post process by combination */
    for(int m = 0; m < M; m++){
        symAfterPP(m) = ComplexD(0,0);
        for(int nr = 0; nr < Nr; ++nr){
            symAfterPP(m) += conj(h_sum[nr]) * symAfterFC(nr, m);
        }
        symAfterPP(m) /= h_conj_sum;
    }

    for(int m = 0; m < M; ++m){
        decode(m) = symAfterPP(m).real() > 0 ? 1:0;
        BER_TOTAL += (decode(m) != source(m));
    }

    #ifdef DebugMode
        cout<<".................AWGN................."<<endl;
        cout<<tmp<<endl; 
        
        cout<<".................symAfterFC................."<<endl;
        cout<<symAfterFC<<endl; 

        cout<<".................h_conj_sum................."<<endl;
        cout<<h_conj_sum<<endl; 
        
        cout<<".................h_sum................."<<endl;
        for(int nr = 0; nr < Nr; ++nr){
            cout<<h_sum[nr]<<" "; 
        }
        cout<<endl;
        cout<<".................symAfterPP................."<<endl;
        cout<<symAfterPP<<endl; 

        cout<<".................source................."<<endl;
        cout<<source<<endl; 
        cout<<".................decode................."<<endl;
        cout<<decode<<endl; 
    #endif
}

