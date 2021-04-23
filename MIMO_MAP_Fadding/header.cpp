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
CSIMatrix H;                                 /* channel parameters , U*Nt*Nr */
ModuMatrix Constell[Mpoint];
SymAfterFCMatrix ConstellFixed[Mpoint];
SourceMatrix Decode;  


void NormalIO(){

    cout<<"MIMO-MAP"<<endl;
    outfile.open("MIMO-MAP.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"MIMO-MAP"<<endl;
    
    cout<<"Nt: "<<Nt<<setw(20)<<"Nr: "<<Nr<<setw(20)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(20)<<"Nr: "<<Nr<<setw(20)<<"NLoop: "<<NLoop<<endl;
    cout<<"power: "<<power<<endl;   
    outfile<<"power: "<<power<<endl; 


    cout<<"MAP"<<endl;
    outfile<<"MAP"<<endl;

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;

    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;
}


void InitMapMatrix(){

    /* initialize the Masterconstellation by bitset */
    for(int i = 0; i < Mpoint; i++){
        bitset<Nt> bit(i);
        for(int nt = 0; nt < Nt; ++nt){
            Constell[i](nt) = ComplexD(static_cast<double>(2 * bit[nt] - 1), 0);
        } 
    }  
}

/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);
    double EsN0 = 1 * ebN0;
    double N0 = power * Nt / EsN0;
    // double N0 = power / (EsN0*Nt);
    N_Var = N0/2;
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
        modu(i) = sqrt(1.0/static_cast<double>(Nt)) * ComplexD(2*source(i)-1, 0);
    }

    #ifdef DebugMode
        cout<<".................bpsk signal................."<<endl;
        cout<<modu<<endl;
        sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            sum += pow(abs(modu(nt)), 2);
        }
        cout << "weight of modu: "<< sum << endl;
    #endif
}


void FadingChannel(CSIMatrix& h){

    /* channel parameters */
    for(int nt = 0; nt < Nt; ++nt){
        for(int nr = 0; nr < Nr; ++nr){
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
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(ModuMatrix& modu, CSIMatrix& h, ModuMatrix* constell,
              SourceMatrix& source,  SourceMatrix& decode){

    /*** Fadding channel ***/
    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; nr++){
        tmp(nr) = AWGN(N_Var);
    }

    SymAfterFCMatrix symAfterFC;
    symAfterFC =sqrt(static_cast<double>(Nt))*(h * modu + tmp);//

    /*** adjust the master-constellation ***/
    for(int point = 0; point < Mpoint; ++point){
        ConstellFixed[point] = h * constell[point];
    }

    #ifdef DebugMode
        cout<<"...................MasterConstell..................."<<endl;
        for(int i = 0; i < Mpoint; i++){    
            cout<<ConstellFixed[i]<<endl;
        }
    #endif  

    /* MAP */
    double minDist = INT16_MAX;
    int minDistPoint = 0;
    double tmpMinDist = 0;

    for(int point = 0; point < Mpoint; ++point){
        tmpMinDist = 0;
        SymAfterFCMatrix tmpMinus = symAfterFC - ConstellFixed[point];

        for(int nr = 0; nr < Nr; ++nr){
            tmpMinDist += pow(abs(tmpMinus(nr)), 2);
        }
        if(tmpMinDist < minDist){
            minDist = tmpMinDist;
            minDistPoint = point;
        }
    }
    bitset<Nt> bit(minDistPoint);
    for(int nt = 0; nt < Nt; ++nt){
        decode(nt) = bit[nt];
    }

    for(int nt = 0; nt < Nt; ++nt){
        if(decode(nt) != source(nt)) BER_TOTAL++;
    }

    #ifdef DebugMode
        cout<<".................AWGN................."<<endl;
        cout<<tmp<<endl; 
        
        cout<<".................bit................."<<endl;
        for(int nt = 0; nt < Nt; ++nt){
            cout<<bit[nt]<<" ";
        }
        cout<<endl;

        cout<<".................symAfterFC................."<<endl;
        cout<<symAfterFC<<endl; 

        cout<<".................minDistPoint................."<<endl;
        cout<<minDistPoint<<endl; 

        cout<<".................decode................."<<endl;
        cout<<decode<<endl; 
    #endif
}

