#include <math.h>
#include <stdio.h>
__declspec(dllexport)
void TransformacionPtoT(double Pa, double Pb, double* Ta, double* Tb) {
    *Ta = 2 * Pa / (1 + Pa - Pb);
    *Tb = 2 * Pb / (1 + Pb - Pa);
}

__declspec(dllexport)
void TransformacionTtoP(double Ta, double Tb, double* Pa, double* Pb) {
    double denom = 2 - Ta - Tb;
    *Pa = (Ta * (1 - Tb)) / denom;
    *Pb = (Tb * (1 - Ta)) / denom;
}

__declspec(dllexport)
double dPa(double Ta, double Tb, double c, double s, double n) {
    return s * (1 - Ta) * (c * (Ta * Ta + (1 - Ta) * (1 - Tb)) + (1 - c) * n)
         - (1 - s) * Ta * (c * ((Ta + Tb) * (1 - Ta)) + (1 - c) * (1 - n));
}

__declspec(dllexport)
double dPb(double Ta, double Tb, double c, double s, double n) {
    return s * (1 - Tb) * (c * (Tb * Tb + (1 - Tb) * (1 - Ta)) + (1 - c) * (1 - n))
         - (1 - s) * Tb * (c * ((Tb + Ta) * (1 - Tb)) + (1 - c) * n);
}

__declspec(dllexport)
void IntegracionEcBias(double c, double sa, double sb, double n, double* Ta_out, double* Tb_out) {
    FILE *f= fopen("Resultados/Numeric.dat", "w");
    int index=0;
    double h = 0.0001;
    double Ta = *Ta_out, Tb = *Tb_out;
    double TaPrev = 0, TbPrev = 0;
    double Pa, Pb, dPa1, dPb1, dPa2, dPb2, Pa_pred, Pb_pred, Ta_pred, Tb_pred;
    TransformacionTtoP(Ta, Tb, &Pa, &Pb);
    while (  sqrt(((Ta - TaPrev)*(Ta - TaPrev) + (Tb - TbPrev)*(Tb - TbPrev))) > 1e-30 && (Ta)<1.0 && (Tb)<1.0 && Ta >= 0.0 && Tb >=0.0) {
        fprintf(f, "%d\t%lf\t%lf\n",index++, Ta, Tb);
        TaPrev = Ta;
        TbPrev = Tb;
        dPa1 = dPa(Ta, Tb, c, sa, n);
        dPb1 = dPb(Ta, Tb, c, sb, n);
        Pa_pred = Pa + h * dPa1;
        Pb_pred = Pb + h * dPb1;

        TransformacionPtoT(Pa_pred, Pb_pred, &Ta_pred, &Tb_pred);


        dPa2 = dPa(Ta_pred, Tb_pred, c, sa, n);
        dPb2 = dPb(Ta_pred, Tb_pred, c, sb, n);

        Pa += h * (dPa1 + dPa2) / 2.0;
        Pb += h * (dPb1 + dPb2) / 2.0;

        TransformacionPtoT(Pa, Pb, &Ta, &Tb);
    }
    *Ta_out = Ta;
    *Tb_out = Tb;
    fclose(f);
}

__declspec(dllexport)
double dQ(double Ta, double Qa, double Qb, double beta, double c){
    return Ta*1.0/(1+exp(-beta*Qa)) - (1-Ta)*1.0/(1+exp(-beta*Qb)) - Qa-c;
}

__declspec(dllexport)
void QLearningdynamics(double Ta, double Tb, double *Qa, double *Qb, double beta, double c){
    double h=0.0001, tolerance = 1e-30;
    double QaPrev = 10, QbPrev=10, QaPred, QbPred, dQ1, dQ2;
    while(sqrt((*Qa - QaPrev)*(*Qa-QaPrev) + (*Qb - QbPrev)*(*Qb-QbPrev)) > tolerance){
        QaPrev = *Qa;
        QbPrev = *Qb;
        dQ1 = dQ(Ta, *Qa, *Qb, beta, c);
        dQ2 = dQ(Tb, *Qb, *Qa, beta, c);
        QaPred = *Qa + h*dQ1;
        QbPred = *Qb + h*dQ2;
        *Qa += h*(dQ1 + dQ(Ta, QaPred, QbPred, beta, c))/2.0;
        *Qb += h*(dQ2 + dQ(Tb, QbPred, QaPred, beta, c))/2.0;
    }
}


int main(){
    double Ta = 0.2, Tb=0.8, n = 0.2, c=0.9, sa=0.755, sb=0.65;
    IntegracionEcBias(c, sa, sb, n, &Ta, &Tb);
}
