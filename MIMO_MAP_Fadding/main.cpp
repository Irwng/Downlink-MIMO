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
    for(int EbN0dB = minEbN0dB; EbN0dB <= maxEbN0dB; EbN0dB = EbN0dB + step){
        ChannelInitialize(EbN0dB);
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            FadingChannel(H);
            Receiver(Modu, H, Constell, Source, Decode);
        }
        BER = static_cast<double>(BER_TOTAL/(NLoop*Nt));
        cout<<EbN0dB<<setw(20)<<BER<<endl;
        outfile<<EbN0dB<<setw(20)<<BER<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();
    return 0;
}
