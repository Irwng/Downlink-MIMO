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
    
    cout<<"MISO-STBC"<<endl;
    outfile.open("MISO-STBC.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"MISO-STBC"<<endl;
    
    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"Time slots: "<<M<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"Time slots: "<<M<<setw(10)<<"NLoop: "<<NLoop<<endl; 
    
    cout<<"Alamouti-STBC"<<endl;
    outfile<<"Alamouti-STBC"<<endl;

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;

    cout<<"snrdB"<<setw(15)<<"BER"<<endl;
    outfile<<"snrdB"<<setw(15)<<"BER"<<endl;
}


/* key point */
void ChannelInitialize(int snrdB){

    double snr = pow(10, (double)snrdB/10);
    N_Var = Nt *power / snr;
    BER_TOTAL = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int i = 0; i < Nt; i++){
        source(i)= rand()%2;
    }

    #ifdef DebugMode
        cout<<".................source................."<<endl;
        cout<<source<<endl;
    #endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int i = 0; i < Nt; i++){
        modu(i) = ComplexD(2*source(i)-1, 0);
    }
    #ifdef DebugMode
        cout<<".................bpsk signal................."<<endl;
        cout<<modu<<endl;
    #endif
}


void Alamouti(ModuMatrix& modu, STMatrix& st){

    st(0, 1) = -conj(modu(1)); st(0, 0) = modu(0);
    st(1, 1) =  conj(modu(0)); st(1, 0) = modu(1);

    #ifdef DebugMode
        cout<<".................Alamouti................."<<endl;
        cout<<st<<endl;
    #endif
}


void FadingChannel(CSIMatrix& h){

    /* channel parameters */
    for(int nt = 0; nt < Nt; ++nt){
        double GG = sqrt(-2*log(randN() + 0.000000001));
        double B = randN();
        h(nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                         sqrt(0.5) * GG * sin(2*PI*B));
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
    symAfterFC = h * symAfterST + tmp;

    double sum = pow(abs(h(0)), 2) + pow(abs(h(1)), 2);
    /* post process by combination */
    symAfterPP(0) = ( conj(h(0)) * symAfterFC(0) + h(1) * conj(symAfterFC(1)) )/sum;
    symAfterPP(1) = ( conj(h(1)) * (symAfterFC(0)) - h(0) * conj(symAfterFC(1)) )/sum;

    for(int j = 0; j < M; ++j){
        decode(j) = symAfterPP(j).real() > 0 ? 1:0;
        BER_TOTAL+=(decode(j) != source(j));
    }

    #ifdef DebugMode
        cout<<".................AWGN................."<<endl;
        cout<<tmp<<endl; 
        
        cout<<".................symAfterFC................."<<endl;
        cout<<symAfterFC<<endl; 

        cout<<".................symAfterPP................."<<endl;
        cout<<symAfterPP<<endl; 
        
        cout<<".................decode................."<<endl;
        cout<<decode<<endl; 
    #endif
}

