/************************************************************
Author: Wangyi  Version: C++11  Date: 2021/3/17
Theme: Downlink MIMO Beamforming  
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

    if(argc < 2 || Nr!=Nt){
        cout<<"1st argument: Preprocess Algorithms"<<endl;
        cout<<"1: ZF"<<endl;
        cout<<"2: MMSE"<<endl;
        cout<<"3: OSIC-ZF-SINR"<<endl;
        cout<<"4: OSIC-MMSE-SINR"<<endl;
        cout<<"3: OSIC-ZF-SNR"<<endl;

        if(Nr!=Nt) cout<<"Dimension Error: U*Nr!=Nt"<<endl;
        return 0;
    }

    double NUM = NLoop*Nr*BitperSymbol;
    #ifndef DebugMode
        Initialize(argv);
    #endif
    for(int EbN0dB = minEbN0dB; EbN0dB <=  maxEbN0dB; EbN0dB=EbN0dB+step){
        ChannelInitialize(EbN0dB);
        for(int loop = 0; loop < NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            if(*argv[1] == '1'|| *argv[1] == '2')
                Receiver(Modu, SymAfterFC, H, V, Source, Decode, argv);
            else
                Receiver_OSIC(Modu, SymAfterFC, H, Source, Decode, argv);
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
