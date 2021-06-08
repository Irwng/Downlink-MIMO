/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/3/16
Theme: MISO-STBC-FaddingChannel
Channel: fading channel
Moduulaiton: BPSK
Antenna: Nt=2,Nr=1
Decoding: MAP
Time slot(M): 2
***********************************************************/
#include "header.h"

fstream outfile;

int main(){

    srand((unsigned)time(NULL));
    time_t start = time(NULL);

    NormalIO();

    for(int snrdB = MinSNRdB; snrdB <= MaxSNRdB; snrdB = snrdB + Step){
        ChannelInitialize(snrdB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            Alamouti(Modu, SymAfterST);
            FadingChannel(H);
            Receiver(SymAfterST, SymAfterFC, H, SymAfterPP, Source, Decode);
        }

        BER = static_cast<double>(BER_TOTAL/(NLoop*Nt));
        cout<<snrdB<<setw(20)<<BER<<endl;
        outfile<<snrdB<<setw(20)<<BER<<endl;
        cout<<"time(s): "<<time(NULL) - start<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();
    return 0;
}
