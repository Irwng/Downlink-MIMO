/************************************************************
Author: Wangyi  Version: C++11  Date: 2021/05/03
Theme: Downlink MU-MIMO Non-orthogonal Beamforming  
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6
Antenna: Nt= U= 2, 1 receiver antenna per user
Beamforming: Non-orthogonal like power-domain NOMA
MAP Receiver
***********************************************************/
#include "header.h"

int main(int argc, char* argv[]){
   
    srand((unsigned)time(NULL));
    time_t start = time(NULL);

    double NUM = NLoop*Nt*U*BitperSymbol;
    Initialize(MasterConstell);

    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Beamforming(Modu, H, W, SymAfterBF, argv);
            Receiver(SymAfterBF, SymAfterFC, H, SymAfterPP, W, Source, Decode, MasterConstell);
            
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
