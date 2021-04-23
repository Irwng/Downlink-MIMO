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
        cout<<"4: DPC"<<endl;

        if(U*Nr!=Nt) cout<<"Dimension Error: U*Nr!=Nt"<<endl;
        return 0;
    }

    double NUM = NLoop*Nt*BitperSymbol;
    Initialize(argv);

    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        alphasum = 0;
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Beamforming(Modu, H, W, SymAfterBF, argv);
            Receiver(SymAfterBF, SymAfterFC, H, Source, Decode, argv);
            
            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }
        BER = static_cast<double>(BER_TOTAL/NUM);
        cout<<snrdB<<setw(20)<<BER<<endl;
        outfile<<snrdB<<setw(20)<<BER<<endl;
    }
    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();
    return 0;
}
