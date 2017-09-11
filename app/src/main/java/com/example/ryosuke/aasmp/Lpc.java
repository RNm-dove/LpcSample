package com.example.ryosuke.aasmp;

import org.apache.commons.math3.complex.*;

/**
 * Created by ryosuke on 17/08/26.
 */

public class Lpc {

    public Lpc() {

    }

    public double[] lpc(double[] r, int order, double df){
        return levionsonDurbin(autocorr(r, order+1), order, df);
    }

    protected class Formant{
        double first;
        double second;
     }

    /**
     *
     * @param x: 信号
     * @param nlags: 自己相関関数のサイズ（lag=0からnlags-1まで)
     *             引数がなければ(lag=0からlen(x)-1まですべて)
     * @return
     */
    public double[] autocorr(double[] x,int nlags){
        int N = x.length;
        if(nlags==0){
            nlags = N;
        }
        double r[] = new double[nlags];
        for(int lag=0; lag <nlags; lag++){
            r[lag] = 0.0;
            for(int n=0; n<(N-lag); n++){
                r[lag] += x[n] * x[n + lag];
            }
        }
        return r;
    }

    public double[] autocorr(double[] x){
        double[] r = autocorr(x, 0);
        return r;
    }

    /**
     * ｋ次のLPC係数からk+1次のLPC係数を再帰的に計算してLPC係数を求める
     * @param r
     * @param lpcOrder
     * @return
     */
    public double[] levionsonDurbin(double[] r,int lpcOrder, double df){
        double a[] = new double[lpcOrder + 1];   //a[0]は１で固定のためlpcOrder個の係数を得るためには+1が必要
        double e[] = new double[lpcOrder + 1];

        // k＝1 のとき
        a[0] = e[0] = 1.0;
        a[1] = - r[1]/ r[0];
        e[1] = r[0] + r[1] * a[1];


        //kの場合からk+1の場合を再帰的に求める
        for(int k=1;k<lpcOrder;k++){
            //lambdaを更新
            double lam = 0.0;
            for(int j=0;j<k+1;j++){
                lam -= a[j] * r[k+1-j];
            }
            lam /= e[k];

            // aを更新
            //UとVからaを更新
            double[] U = new double[k+2];
            double[] V = new double[k+2];
            U[0] = 1.0;
            V[0] = 0.0;
            for(int i=1;i<k+1;i++){
                U[i] = a[i];
                V[k + 1 -i] = a[i];
            }
            U[k + 1] = 0.0;
            V[k + 1] = 1.0;


            for(int i=0; i< k + 2; i++){
                a[i] = U[i] + lam*V[i];
            }

            //eを更新
            e[k+1] = e[k] * (1.0 - lam * lam);

        }

        return freqz(e, a, df, r.length);

    }

    protected double[] freqz(double[] e, double[] a, double df, int N){
        double[] H = new double[N];
        for(int n=0; n < N; ++n){

            Complex w = new Complex(0.0, -2.0 * Math.PI * n/N );
            Complex z = w.exp();
            Complex numerator = new Complex(0.0, 0.0);
            Complex denominator = new Complex(0.0, 0.0);

            for(int i = 0; i< e.length ; ++i){
                numerator = numerator.add(z.pow(i).multiply(e[e.length - 1 -i]));
            }
            for (int i =0; i < a.length; ++i){
                denominator = denominator.add(z.pow(i).multiply(a[a.length - 1 -i]));
            }
            H[n] = numerator.divide(denominator).abs();

        }

        return H;

    }

    public Formant formant(double[] r, double df){
        Formant result = new Formant();
        result.first = 0.0;
        result.second = 0.0;
        boolean is_find_first = false;
        for(int i = 1; i< r.length -1; ++i){
            if(r[i] > r[i-1] && r[i]>r[i+1]){
                if(!is_find_first){
                    result.first = df * i;
                    is_find_first = true;
                } else {
                    result.second = df * i;
                    break;
                }
            }
        }

        return result;

    }

    public double[] normalize(double[] r){
        double max = r[0];
        double min = r[1];
        for(int i = 1; i < r.length; i++){
            double v = r[i];
            if(v > max){
                max = v;
            }
            if(v < min){
                min = v;
            }
        }
        double factor = Math.max(Math.abs(max),Math.abs(min));
        double[] result = new double[r.length];
        for(int i=0; i<r.length; ++i){
            result[i] = r[i]/factor;
        }
        return result;
    }

    public double[] preEmphasis(short[] r){
        double[] result = new double[r.length];
        result[0] = (double)r[0];
        for(int i=1; i< r.length;i++){
            result[i] = (double)r[i] - 0.98*r[i-1];
        }
        return result;
    }

    public double[] hamming(double[] r){
        int N = r.length;
        double[] result = new double[N];
        for(int i=1; i<N -1; ++i){
            double h = 0.54 - 0.46 * Math.cos(2* Math.PI * i / (N -1));
            result[i] = r[i] * h;
        }
        result[0] = result[N-1] = 0;
        return result;
    }

    public short[] toShort(double[] r){
        short[] result = new short[r.length];
        for(int i=0 ; i<r.length; i++){
            result[i] = (short)(32767*r[i]);
        }
        return result;
    }

    public double volume(double[] r){
        double v = 0.0;
        for(double x: r){
            v += x*x;
        }
        v /= r.length;
        return v;
    }

    public String vowel(double f1, double f2)
    {
        if (f1 > 600 && f1 < 1400 && f2 > 900  && f2 < 2000) return "あ";
        if (f1 > 100 && f1 < 410  && f2 > 1900 && f2 < 3500) return "い";
        if (f1 > 100 && f1 < 700  && f2 > 1100 && f2 < 2000) return "う";
        if (f1 > 400 && f1 < 800  && f2 > 1700 && f2 < 3000) return "え";
        if (f1 > 300 && f1 < 900  && f2 > 500  && f2 < 1300) return "お";
        return "-";
    }


}
