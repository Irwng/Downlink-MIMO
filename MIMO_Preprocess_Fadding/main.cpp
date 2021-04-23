/************************************************************
Author: Wangyi  Version: C++11  Date: 2021/3/17
Theme: Downlink MIMO Beamforming  
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6
Antenna: Nt= U= 2, 2 receiver antenna per user
Beamforming: ZF/MMSE
Note: 2021/4/02- realize the OSIC algorithm getting better performance
      2021/4/17- modify the OSIC algorithm saving half time
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
        cout<<"5: OSIC-ZF-SNR"<<endl;
        cout<<"6: OSIC-MMSE-SNR"<<endl;

        if(Nr!=Nt) cout<<"Dimension Error: U*Nr!=Nt"<<endl;
        return 0;
    }

    double NUM = NLoop*Nt*BitperSymbol;
    #ifndef DebugMode
        Initialize(argv);
    #endif
    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            if(*argv[1] == '1'|| *argv[1] == '2'){
                Receiver(Modu, SymAfterFC, H, V, Source, Decode, argv);
            }
            else{
                Receiver_OSIC(Modu, SymAfterFC, H, Source, Decode, argv);
                // Receiver_OSIC_V2(Modu, SymAfterFC, H, V, Source, Decode, argv);
            }
           
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
