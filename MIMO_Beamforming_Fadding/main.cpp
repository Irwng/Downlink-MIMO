/************************************************************
Author: Wangyi  Version: C++11  Date: 2021/3/17
Theme: Downlink MU-MIMO Beamforming  
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6
Antenna: Nt= U= 2, 2 receiver antenna per user
Beamforming: ZF/MMSE
***********************************************************/
#include "header.h"

int main(int argc, char* argv[]){
   
    srand((unsigned)time(NULL));
    time_t start = time(NULL);

    if(argc < 2 || U*Nr!=Nt){
        cout<<"1st argument: Beamforming Algorithms"<<endl;
        cout<<"1: ZF"<<endl;
        cout<<"2: MMSE"<<endl;
        cout<<"3: SVD"<<endl;
        cout<<"4: BD"<<endl;

        if(U*Nr!=Nt) cout<<"Dimension Error: U*Nr!=Nt"<<endl;
        return 0;
    }

    double NUM = NLoop*Nr*U*BitperSymbol;
    Initialize(argv);

    for(int EbN0dB = minEbN0dB; EbN0dB <= maxEbN0dB; EbN0dB=EbN0dB+2){
        ChannelInitialize(EbN0dB);
        alphasum = 0;
        for(int loop = 0; loop < NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Beamforming(Modu, H, W, SymAfterBF, argv);
            Receiver(SymAfterBF, SymAfterFC, H, Source, Decode, argv);
        }
        BER = static_cast<double>(BER_TOTAL/NUM);
        cout<<EbN0dB<<setw(20)<<BER<<endl;
        outfile<<EbN0dB<<setw(20)<<BER<<endl;
    }
    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();
    return 0;
}
