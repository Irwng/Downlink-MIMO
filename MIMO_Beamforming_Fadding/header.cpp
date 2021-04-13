/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;                                     /* Bits Error Rate */
double alpha = 0;
double alphasum = 0;
fstream outfile;

SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
CSIMatrix H;                                        
CSIMatrix svdU;                                     
CSIMatrix svdV;                                     
ModuMatrix vectorS;
CSIMatrix WeightedIdentityMatrix;                   /* MMSE assistance matrix */
BFMatrix W;                                         /* beamforming matrix */
SymAfterFCMatrix SymAfterFC;                        /* receiving signals after postprocessing */
SymAfterBFMatrix SymAfterBF;                        /* receiving signals after postprocessing */
DecodeMatrix Decode;  


void Initialize(char* argv[]){

    outfile.open("Downlink MU-MIMO Beamforming.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    cout<<"Downlink MU-MIMO Beamforming"<<endl;
    outfile<<"Downlink MU-MIMO Beamforming"<<endl;

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
            cout<<"SVD"<<endl;
            outfile<<"SVD"<<endl;
            break;

        case '4':
            cout<<"BD"<<endl;
            outfile<<"BD"<<endl;
            break;

        default:
            cout<<"Invaild input!"<<endl;
            break;
    }

    cout<<"FaddingChannel"<<endl;
    outfile<<"FaddingChannel"<<endl;
    
    cout<<"Order of modulation: "<<Mod<<endl;   
    outfile<<"Order of modulation: "<<Mod<<endl; 
    
    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl; 
   
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;
}


void ChannelInitialize(int ebN0dB){

    /* the total power of transmitter is fixed */
    double ebN0 = pow(10,(double)ebN0dB/10);
    double snr = BitperSymbol * ebN0;
    double N0 = power * Nt /(snr);
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
    for(int i = 0; i < LenBit; i++){
        source(i)= rand()%2;
    }
    #ifdef DebugMode
        cout<<"source: "<<endl<<source<<endl;
    #endif
}


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    #ifdef BPSKMode
        /* BPSK: 0->-1, 1->1 */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD(2 * source(nt) - 1, 0);
        }
    #else
        /* QPSK: 00->a+aj, 01->-a+aj, 11-> -a-aj, 10->a-aj */
        /* QPSK: forward bits->1st antenna, next bits->2st antenna */
        for(int nt = 0; nt < Nt; ++nt){
            modu(nt) = ComplexD((1 - 2 * source(BitperSymbol * nt + 1)) * sqrt(0.5), (1 - 2 * source(BitperSymbol * nt)) * sqrt(0.5));
        }
    #endif

    #ifdef DebugMode
        cout<<".................Alamouti................."<<endl;
        cout<<modu<<endl;
        double sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            sum += pow(abs(modu(nt)), 2);
        }
        cout << "weight of modu: "<< sum << endl;
    #endif


    #ifdef DebugMode
       cout<<"bpsk modulation: "<<endl<<modu<<endl;
    #endif  
}


void FadingChannel(CSIMatrix& h){

    /* channel parameters */
    for(int u = 0; u < U*Nr; ++u){
        for(int nt = 0; nt < Nt; ++nt){
            double GG = sqrt(-2*log(randN() + 0.000000001));
            double B = randN();
            h(u, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
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


void Beamforming(ModuMatrix& modu, CSIMatrix& h,
                 BFMatrix& w, SymAfterBFMatrix& symAfterBF,
                 char* argv[]){

    CSIMatrix Denomiator_H;
    CSIMatrix ConjTrans_H;
    CSIMatrix Base_H;
    double sum = 0;
    
    switch (*argv[1]){
        /* ZF */
        case '1':
            ConjTrans_H = h.conjugate().transpose();
            Denomiator_H = h * ConjTrans_H; 
            Base_H = ConjTrans_H * (Denomiator_H.inverse());
            
            w = Base_H;
            symAfterBF = w * modu;

            sum = 0;
            for(int nt = 0; nt < Nt; ++nt){
                sum += pow(abs(symAfterBF(nt)), 2);
            }
            symAfterBF = sqrt(power/sum) * symAfterBF;

            break;

        /* MMSE */
        case '2':
            ConjTrans_H = h.conjugate().transpose();
            Denomiator_H = h * ConjTrans_H + WeightedIdentityMatrix; 
            Base_H = ConjTrans_H * (Denomiator_H.inverse());

            w = Base_H;
            symAfterBF = w * modu;

            sum = 0;
            for(int nt = 0; nt < Nt; ++nt){
                sum += pow(abs(symAfterBF(nt)), 2);
            }
            symAfterBF = sqrt(power/sum) * symAfterBF;

            break;

        /* SVD */
        case '3':
            SVD(h, svdU, svdV, vectorS);
            symAfterBF =  sqrt(0.5) * svdV * modu;
            break;

        default:
            break;
    }

    #ifdef DebugMode
        cout<<"w"<<endl<<w<<endl;
        sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            for(int nr = 0; nr < Nt; ++nr){
                sum += pow(abs(w(nt, nr)), 2);
            }
        }
        cout << "weight of w: "<< sum << endl;
        
        cout<<"symAfterBF"<<endl<<symAfterBF<<endl;
        sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            sum += pow(abs(symAfterBF(nt)), 2);
        }
        cout<<"weight of symAfterBF: "<<sum<<endl;
    #endif
}


void Receiver(SymAfterBFMatrix& symAfterBF, SymAfterFCMatrix& symAfterFC, 
              CSIMatrix& h, SourceMatrix& source,
              DecodeMatrix& decode, char* argv[]){

    /* AWGN */
    SymAfterFCMatrix tmp;
    for(int nr = 0; nr < U*Nr; ++nr){
        tmp(nr) = AWGN(N_Var);
    }

    switch (*argv[1]){
        case '1':
        case '2':
            symAfterFC = h * symAfterBF + tmp;
            break;

        case '3':
            symAfterFC = (svdU.conjugate().transpose())* h * symAfterBF + tmp;
            break;

        default:
            break;
    }

    /* check the performance of the post processing */
    #ifdef BPSKMode
        for(int nr = 0; nr < U*Nr; ++nr){
            decode(nr) = static_cast<int>(symAfterFC(nr).real() > 0);
            BER_TOTAL += static_cast<double>(decode(nr) != source(nr));
        }
    #else
        for(int nr = 0; nr < U*Nr; ++nr){
            decode(nr * BitperSymbol) = static_cast<int>(symAfterFC(nr).imag() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol)!=source(nr * BitperSymbol));
            decode(nr * BitperSymbol + 1) = static_cast<int>(symAfterFC(nr).real() < 0);
            BER_TOTAL += static_cast<double>(decode(nr * BitperSymbol + 1)!=source(nr * BitperSymbol + 1));
        }
    #endif

}


void SVD(CSIMatrix& h, CSIMatrix& u, CSIMatrix& v, ModuMatrix& vectorS){

    JacobiSVD<CSIMatrix> svd(h, ComputeFullU | ComputeFullV);
    u = svd.matrixU();
    v = svd.matrixV();
    vectorS = svd.singularValues();

    #ifdef DebugMode
        CSIMatrix tmpU = u * (u.conjugate().transpose());
        CSIMatrix tmpV = v * (v.conjugate().transpose());
        // CSIMatrix S =  u.inverse() * h * (v.conjugate().transpose().inverse()); // S = U^-1 * A * VT * -1
        cout << "JacobiSVD---------------------------------------------"<< std::endl;
        cout << "H"<< endl << h << endl;
        cout << "U*U^H"<< endl << tmpU << endl;
        cout << "V*V^H"<< endl << tmpV << endl;
        // cout << "vectorS"<< endl << vectorS << endl;

        double sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            for(int nr = 0; nr < Nt; ++nr){
                sum += pow(abs(tmpU(nt, nr)), 2);
            }
        }
        cout << "weight of U*U^H: "<< sum << endl;

        sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            for(int nr = 0; nr < Nt; ++nr){
                sum += pow(abs(tmpV(nt, nr)), 2);
            }
        }
        cout << "weight of V*V^H: "<< sum << endl;
        
    #endif
}