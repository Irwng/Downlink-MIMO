/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/3/16
Theme: MIMO-MAP-FaddingChannel
Channel: fading channel
Moduulaiton: BPSK
Decoding: MAP
Note:the number of antennas can be adjusted 
***********************************************************/
#include "header.h"

fstream outfile;

int main(){

    srand((unsigned)time(NULL));
    time_t start = time(NULL);

    NormalIO();

    InitMapMatrix();
    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Receiver(SymAfterFC, Modu, H, Constell, Source, Decode);
            
            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }
        BER = static_cast<double>(BER_TOTAL/(NLoop*Nt));
        cout<<snrdB<<setw(20)<<BER<<endl;
        outfile<<snrdB<<setw(20)<<BER<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();
    return 0;
}
