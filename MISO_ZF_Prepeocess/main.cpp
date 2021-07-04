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

    double NUM = NLoop*Nt*BitperSymbol;
    Initialize();
    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Receiver(Modu, SymAfterFC, H, V, Source, Decode, SymAfterPP);
           
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
