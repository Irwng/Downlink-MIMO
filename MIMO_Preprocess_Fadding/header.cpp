/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;                                     /* Bits Error Rate */
fstream outfile;

SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
CSIMatrix H;                                        
CSIMatrix WeightedIdentityMatrix;                   /* MMSE assistance matrix */
CSIMatrix V;                                        /* preprocess matrix */
SymAfterFCMatrix SymAfterFC;                        /* receiving signals */ 
SymAfterFCMatrix SymAfterPP;                        /* receiving signals after postprocessing */
DecodeMatrix Decode;  


void Initialize(char* argv[]){

    outfile.open("Downlink MIMO Preprocess.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    cout<<"Downlink MIMO Preprocess"<<endl;
    outfile<<"Downlink MIMO Preprocess"<<endl;

    switch (*argv[1]){
    
        case '1':
            cout<<"ZF"<<endl;
            outfile<<"ZF"<<endl;
            break;

        case '2':
            cout<<"MMSE"<<endl;
            outfile<<"MMSE"<<endl;
            break;
            
        case '3':
            cout<<"OSIC-ZF-SINR"<<endl;
            outfile<<"OSIC-ZF-SINR"<<endl;
            break;

        case '4':
            cout<<"OSIC-MMSE-SINR"<<endl;
            outfile<<"OSIC-MMSE-SINR"<<endl;
            break;

        case '5':
            cout<<"OSIC-ZF-SNR"<<endl;
            outfile<<"OSIC-ZF-SNR"<<endl;
            break;

        default:
            cout<<"Invaild input!"<<endl;
            break;
    }

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;

    #ifdef BPSK
        cout<<"BPSK"<<endl;   
        outfile<<"BPSK"<<endl; 
    #endif
    #ifdef QPSK 
        cout<<"QPSK"<<endl;   
        outfile<<"QPSK"<<endl;			            
    #endif
    #ifdef QAM16
        cout<<"QAM16"<<endl;   
        outfile<<"QAM16"<<endl;			            
    #endif

    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"NLoop: "<<NLoop<<endl; 
   
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;
}


void ChannelInitialize(int ebN0dB){

    /* the total power of transmitter is fixed */
    double ebN0 = pow(10,(double)ebN0dB/10);
    double snr = BitperSymbol * ebN0;
    double N0 = Nt / snr;
    N_Var = N0/2;
    
    for(int nt = 0; nt < Nt; ++nt){
        WeightedIdentityMatrix(nt, nt) = ComplexD(N_Var, 0);
    }
    
    BER_TOTAL = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
        cout<<"WeightedIdentityMatrix"<<endl<<WeightedIdentityMatrix<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int nt = 0; nt < Nt; ++nt){
        source(nt)= rand()%2;
    }
    #ifdef DebugMode
        cout<<"source: "<<endl<<source<<endl;
    #endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    #ifdef BPSK
        /* BPSK: 0->-1, 1->1 */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD(2 * source(nt) - 1, 0);
        }
        double sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            sum += pow(abs(modu(nt)), 2);
        }
        modu = sqrt(power/sum) * modu;
    #endif
    #ifdef QPSK
        /* QPSK: 00->a+aj, 01->-a+aj, 11-> -a-aj, 10->a-aj */
        /* QPSK: forward bits->1st antenna, next bits->2st antenna */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD((1 - 2 * source(BitperSymbol * nt + 1)) * sqrt(0.5), (1 - 2 * source(BitperSymbol * nt)) * sqrt(0.5));
        }
    #endif

    #ifdef QAM16
        /* 16QAM: 00->a+aj, 01->-a+aj, 11-> -a-aj, 10->a-aj */
        /* QPSK: forward bits->1st antenna, next bits->2st antenna */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD((1 - 2 * source(BitperSymbol * nt + 1)) * sqrt(0.5), (1 - 2 * source(BitperSymbol * nt)) * sqrt(0.5));
        }
    #endif

    
    #ifdef DebugMode
       cout<<"modulation: "<<endl<<modu<<endl;
    #endif  
}


void FadingChannel(CSIMatrix& h){

    /* channel parameters */
    for(int nr = 0; nr < Nr; ++nr){
        for(int nt = 0; nt < Nt; ++nt){
            double GG = sqrt(-2*log(randN() + 0.000000001));
            double B = randN();
            h(nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                sqrt(0.5) * GG * sin(2*PI*B));
        }
    }
    #ifdef DebugMode
        cout<<"h"<<endl<<h<<endl;
    #endif 
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(ModuMatrix& modu, SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h, CSIMatrix& v, SourceMatrix& source, 
              DecodeMatrix& decode, char* argv[]){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; ++nr){
        tmp(nr) = AWGN(N_Var);
    }
    symAfterFC = h * modu + tmp;

    CSIMatrix Denomiator_H;
    CSIMatrix ConjTrans_H;

    ConjTrans_H = h.conjugate().transpose();
    
    switch (*argv[1]){
        case '1':
            Denomiator_H = ConjTrans_H * h;
            break;
        
        case '2':
            Denomiator_H = ConjTrans_H * h + WeightedIdentityMatrix;
            break;

        default:
            break;
    }

    v = Denomiator_H.inverse() * ConjTrans_H;
    symAfterFC = v * symAfterFC;

    /* check the performance of the post processing */
    #ifdef BPSK
        for(int nr = 0; nr < Nr; ++nr){
            decode(nr) = static_cast<int>(symAfterFC(nr).real() > 0);
            BER_TOTAL += static_cast<double>(decode(nr) != source(nr));
        }
    #endif
    #ifdef QPSK
        for(int nr = 0; nr < Nr; ++nr){
            decode(nr * BitperSymbol) = static_cast<int>(symAfterFC(nr).imag() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol)!=source(nr * BitperSymbol));
            decode(nr * BitperSymbol + 1) = static_cast<int>(symAfterFC(nr).real() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol + 1)!=source(nr * BitperSymbol + 1));
        }
    #endif

    #ifdef QAM16
        for(int nr = 0; nr < Nr; ++nr){
            decode(nr * BitperSymbol) = static_cast<int>(symAfterFC(nr).imag() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol)!=source(nr * BitperSymbol));
            decode(nr * BitperSymbol + 1) = static_cast<int>(symAfterFC(nr).real() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol + 1)!=source(nr * BitperSymbol + 1));
        }
    #endif
}

void Receiver_OSIC(ModuMatrix& modu, SymAfterFCMatrix& symAfterFC, 
                   CSIMatrix& h, SourceMatrix& source, 
                   DecodeMatrix& decode, char* argv[]){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < Nr; ++nr){
        tmp(nr) = AWGN(N_Var);
    }
    symAfterFC = h * modu + tmp;//  
    
    bool index_array[Nt]; // index of the decoded signal
    for(auto& a: index_array) a = false;
    for(int nt = 0; nt < Nt; ++nt){
        /* update the channel gain matrix */
        ComplexD* h_buff = new ComplexD[Nt*(Nt-nt)];
        
        int column = 0;
        for(int col = 0; col < Nt; ++col){
            if(index_array[col] == true) continue;
            for(int row = 0; row < Nt; ++row){
                h_buff[column*Nt + row] = h(row, col);
            }        
            column++;
        }
        
        Matrix<ComplexD, Dynamic, Dynamic> h_temp;
        h_temp = Map<Matrix<ComplexD, Dynamic, Dynamic>>(h_buff, Nt, Nt-nt);                
        
        #ifdef DebugMode
            cout<<"nt: "<<nt<<endl;
            cout<<"h_buff: "<<endl;
            for(int k = 0; k < Nt*(Nt-nt); ++k){
                cout<<h_buff[k]<<" ";
            }
            cout<<endl;
            cout<<"h_temp:"<<endl<<h_temp<<endl;
        #endif
        delete[] h_buff;

        /* calculate the h_temp's conjugate transpose matrix */
        Matrix<ComplexD, Dynamic, Dynamic> ConjTrans_H;
        ConjTrans_H = h_temp.conjugate().transpose();

        /* calculate the base matrix */
        Matrix<ComplexD, Dynamic, Dynamic> Denomiator_H;
        Denomiator_H = ConjTrans_H * h_temp;
        switch (*argv[1]){
            case '3':
            case '5':
                break;
            case '4':
                for(int i = 0; i < Nt-nt; ++i) Denomiator_H(i, i) += ComplexD(N_Var, 0);
                break;                    
            default:
                break;
        }
        /* calculate the preprocess matrix */
        Matrix<ComplexD, Dynamic, Dynamic> V_temp;
        V_temp = Denomiator_H.inverse() * ConjTrans_H;

        Matrix<ComplexD, Dynamic, Dynamic> VzfH;
        VzfH = V_temp * h_temp;

        /* calculate the parameters need by SINR */
        double powersum[Nt-nt];
        double noisesum[Nt-nt];
        for(int i = 0; i < Nt-nt; ++i){
            powersum[i] = 0;
            noisesum[i] = 0;
            for(int j = 0; j < Nt-nt; ++j)
                powersum[i] += pow(abs(VzfH(i, j)), 2);

            for(int j = 0; j < Nt; ++j)
                noisesum[i] += pow(abs(V_temp(i, j)), 2);
        }
        
        float SINR[Nt-nt] = {0};
        float denominator = 0.0;
        
        /* calculate the maximum SINR */        
        switch (*argv[1]){
            case '3':
            case '4':
                for(int nr = 0; nr < Nr-nt; ++nr){
                    denominator = powersum[nr] - pow(abs(VzfH(nr,nr)), 2) + N_Var * noisesum[nr];
                    SINR[nr] = pow(abs(VzfH(nr,nr)), 2) / denominator;
                }
                break;
            case '5':
                for(int nr = 0; nr < Nr-nt; ++nr){
                    denominator = N_Var * noisesum[nr];
                    SINR[nr] = 1 / denominator;
                }
                break;                    
            default:
                break;
        }
        #ifdef DebugMode
            cout<<"SINR"<<endl;
            for(auto b: SINR) cout<<b<<" ";
            cout<<endl;
        #endif

        /* find the max SINR subchannel */
        int index_max = 0;
        for(int i = 1; i < Nt-nt; ++i){
            if(SINR[i] > SINR[index_max]) index_max = i;         
        }

        /* find the order in ture line*/
        /* initialize */
        for(int i = 0; i < Nt; ++i){
            if(index_array[i] == false){
                column = i;
                break;    
            }
        }

        for(int i = 0; i < index_max; ++i){
            column++;
            while(column < Nt&& index_array[column] == true){
                column++;
            }
        }

        if(index_array[column] == true) cout<<"error!";
        index_array[column] = true;

        #ifdef DebugMode
            cout<<"index_max: "<<index_max<<endl;
            cout<<"column: "<<column<<endl;
            cout<<"index_array:"<<endl;
            for(auto b: index_array) cout<<b<<" ";
            cout<<endl;
        #endif

        /* decode the signal on the max SINR subchannel in this loop */
        ComplexD x_temp(0,0);
        for(int nr = 0; nr < Nt; ++nr) x_temp += V_temp(index_max, nr) * symAfterFC(nr);
        /* check the performance of the post processing */
        #ifdef BPSK
            decode(column) = static_cast<int>(x_temp.real() > 0);
            x_temp = ComplexD(static_cast<double>(2*decode(column) - 1),0);
        #endif
        #ifdef QPSK
            decode(column * BitperSymbol) = static_cast<int>(x_temp.imag() < 0);
            decode(column * BitperSymbol + 1) = static_cast<int>(x_temp.real() < 0);
            x_temp = ComplexD((1 - 2 * decode(column * BitperSymbol + 1)) * sqrt(0.5), (1 - 2 * decode(column * BitperSymbol)) * sqrt(0.5));

        #endif
        #ifdef QAM16
            decode(column * BitperSymbol) = static_cast<int>(x_temp.imag() < 0);
            decode(column * BitperSymbol + 1) = static_cast<int>(x_temp.real() < 0);
        #endif
        
        /* cancell the interference */
        for(int nr = 0; nr < Nt; ++nr) symAfterFC(nr) -= h_temp(nr, index_max) *sqrt(1.0/Nt) *x_temp;
    }
 
    /* check the performance of the post processing */
    #ifdef BPSK
        for(int nr = 0; nr < Nr; ++nr){
            BER_TOTAL += static_cast<double>(decode(nr) != source(nr));
        }
    #endif
    #ifdef QPSK
        for(int nr = 0; nr < Nr; ++nr){
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol)!=source(nr * BitperSymbol));
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol + 1)!=source(nr * BitperSymbol + 1));
        }
    #endif
    #ifdef QAM16
        for(int nr = 0; nr < Nr; ++nr){
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol)!=source(nr * BitperSymbol));
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol + 1)!=source(nr * BitperSymbol + 1));
        }
    #endif
}
