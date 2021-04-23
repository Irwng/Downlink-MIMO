/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"
double N_Var;
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;                                     /* Bits Error Rate */
double power;
double belta;
double alpha = 0;
double alphasum = 0;
fstream outfile;

SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
CSIMatrix H;                                        
CSIMatrix svdU;                                     
CSIMatrix svdV; 
CSIMatrix qrQ;
CSIMatrix qrR;                                    
CSIMatrix qrA;                                    
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
            cout<<"DPC"<<endl;
            outfile<<"DPC"<<endl;
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
    
    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"U: "<<U<<setw(10)<<"NLoop: "<<NLoop<<endl; 
   
    cout<<"SNRdB"<<setw(15)<<"BER"<<endl;
    outfile<<"SNRdB"<<setw(15)<<"BER"<<endl;
}


void ChannelInitialize(int snrdB){

    /* the total power of transmitter is fixed */
    double snr = pow(10,(double)snrdB/10);
    power = 1.0;
    N_Var = 1.0 / snr;
    
    for(int nt = 0; nt < Nt; ++nt){
        for(int nr = 0; nr < Nt; ++nr){
            WeightedIdentityMatrix(nr, nt) = ComplexD(0, 0);
        }
    }

    for(int nt = 0; nt < Nt; ++nt){
        WeightedIdentityMatrix(nt, nt) = ComplexD(N_Var*Nt, 0);
    }
    
    for(int nt = 0; nt < Nt; ++nt){
        for(int nr = 0; nr < Nt; ++nr){
            qrA(nt, nr) = ComplexD(0, 0);
        }
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
        cout<<modu<<endl;
        double sum = 0;
        for(int nt = 0; nt < Nt; ++nt){
            sum += pow(abs(modu(nt)), 2);
        }
        cout << "weight of modu: "<< sum << endl;
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
    ComplexD GuassN(sqrt(nvar/2.0)* GG * cos(2*PI*B), sqrt(nvar/2.0)* GG * sin(2*PI*B));
    return GuassN;
}


void Beamforming(ModuMatrix& modu, CSIMatrix& h, BFMatrix& w, 
                 SymAfterBFMatrix& symAfterBF, char* argv[]){

    CSIMatrix Denomiator_H;
    CSIMatrix ConjTrans_H;
    CSIMatrix Base_H;
    CSIMatrix Denomiator_W;
    double sum = 0;
    
    switch (*argv[1]){
        /* ZF */
        case '1':
            ConjTrans_H = h.conjugate().transpose();
            Denomiator_H = h * ConjTrans_H; 
            Base_H = ConjTrans_H * (Denomiator_H.inverse());

            Denomiator_W = Base_H * (Base_H.conjugate().transpose());
            sum = 0;
            for(int nt = 0; nt < Nt; ++nt){
                sum += abs(Denomiator_W(nt,nt));
            }
            belta = sqrt(power/sum); 

            w = belta * Base_H;
            symAfterBF = w * modu;

            break;

        /* MMSE */
        case '2':
            ConjTrans_H = h.conjugate().transpose();
            Denomiator_H = h * ConjTrans_H + WeightedIdentityMatrix; 
            Base_H = ConjTrans_H * (Denomiator_H.inverse());

            Denomiator_W = Base_H * (Base_H.conjugate().transpose());
            sum = 0;
            for(int nt = 0; nt < Nt; ++nt){
                sum += abs(Denomiator_W(nt,nt));
            }
            belta = sqrt(power/sum); 
            w = belta * Base_H;
            symAfterBF = w * modu;

            break;

        /* SVD */
        case '3':
            SVD(h, svdU, svdV, vectorS);
            symAfterBF = svdV * modu;

            break;

        /* DPC */
        case '4':
            QR(h, qrQ, qrR);
            for(int nt = 0; nt < Nt; ++nt) qrA(nt, nt) = ComplexD(1, 0) /qrR(nt, nt);
            w = qrQ * qrA;
            
            Base_H = w * ((h * w).inverse());

            Denomiator_W = Base_H * (Base_H.conjugate().transpose());
            sum = 0;
            for(int nt = 0; nt < Nt; ++nt){
                sum += abs(Denomiator_W(nt,nt));
            }
            belta = sqrt(power/sum);
            w = belta * Base_H;
            symAfterBF = w * modu;

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
        case '4':
            // symAfterFC = (1.0/belta)*(h * symAfterBF + tmp);
            symAfterFC = h * symAfterBF + tmp;
            break;

        case '3':
            symAfterFC = (svdU.conjugate().transpose()) * (h * symAfterBF + tmp);//
            for(int nr = 0; nr < Nt; ++nr){
                symAfterFC(nr) /= vectorS(nr);
            }
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

#ifdef DebugMode
    cout<<"symAfterFC"<<endl<<symAfterFC<<endl;
#endif

}


void SVD(CSIMatrix& h, CSIMatrix& u, CSIMatrix& v, ModuMatrix& vectorS){

    JacobiSVD<CSIMatrix> svd(h, ComputeFullU | ComputeFullV);
    u = svd.matrixU();
    v = svd.matrixV();
    vectorS = svd.singularValues();

    #ifdef DebugMode
        CSIMatrix S =  u.inverse() * h * (v.conjugate().transpose().inverse()); // S = U^-1 * H * VT^-1
        cout << "---------------------JacobiSVD------------------------"<< std::endl;
        cout << "H"<< endl << h << endl;
        cout << "S"<< endl << S << endl;
        cout << "U*S*V^H"<< endl << u*S*(v.conjugate().transpose()) << endl;
    #endif
}

void QR(CSIMatrix& h, CSIMatrix& q, CSIMatrix& r){

    HouseholderQR<CSIMatrix> qr;
    qr.compute(h.conjugate().transpose());
    r = qr.matrixQR().triangularView<Eigen::Upper>();
    q = qr.householderQ();

#ifdef DebugMode
    cout << "H = " << endl << h << endl << endl;
    cout << "Q = " << endl << q << endl << endl;
    cout << "R = " << endl << r << endl << endl;
#endif
}