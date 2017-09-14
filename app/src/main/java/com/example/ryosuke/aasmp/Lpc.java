package com.example.ryosuke.aasmp;

import android.util.Log;

import org.apache.commons.math3.complex.*;

import java.util.Arrays;

/**
 * Created by ryosuke on 17/08/26.
 */

public class Lpc {

    public Lpc() {

    }

    public double[] lpc(double[] input, int order, double df){

        int N = input.length;

        //自己相関関数
        double[] r = new double[N];

        for(int l=0; l <= order; l++ ){
            double d = 0.0;
            for(int n=0; n < N-l; n++){
                d += input[n] * input[n + l];
            }
            r[l] = d;
        }

        //levinson-durbin
        double a[] = new double[order + 1];   //a[0]は１で固定のためlpcOrder個の係数を得るためには+1が必要
        double e[] = new double[order + 1];
        Arrays.fill(a, 0);
        Arrays.fill(e, 0);

        // k＝1 のとき
        a[0] = e[0] = 1.0;
        a[1] = - r[1]/ r[0];
        e[1] = r[0] + r[1] * a[1];

        //kの場合からk+1の場合を再帰的に求める
        for(int k=1; k < order; k++){
            //lambdaを更新
            double lam = 0.0;
            for(int j=0; j < k + 1; j++){
                lam -= a[j] * r[k + 1 - j];
            }
            lam /= e[k];

            // aを更新
            //UとVからaを更新
            double[] U = new double[k+2];
            double[] V = new double[k+2];
            U[0] = 1.0;
            V[0] = 0.0;
            for(int i = 1; i < k + 1; i++){
                U[i] = a[i];
                V[k + 1 - i] = a[i];
            }
            U[k + 1] = 0.0;
            V[k + 1] = 1.0;


            for(int i=0; i< k + 2; i++){
                a[i] = U[i] + lam * V[i];
            }

            //eを更新
            e[k+1] = e[k] * (1.0 - lam * lam);

        }

        return  freqz(e, a, df, N);



    }

    protected class Formant{
        double first;
        double second;
     }


    protected double[] freqz(double[] e, double[] a, double df, int N){
        double[] H = new double[N];
        for(int n = 0; n < N  ; n++){

            Complex w = new Complex(0.0, -2.0 * Math.PI * (double)n/N );
            Complex z = w.exp();
            Complex numerator = new Complex(0.0, 0.0);
            Complex denominator = new Complex(0.0, 0.0);

            for(int i = 0; i< e.length ; ++i){
                numerator = numerator.add(z.pow(i).multiply(e[e.length - 1 -i]));
            }
            for (int i =0; i < a.length; ++i){
                denominator = denominator.add(z.pow(i).multiply(a[a.length - 1 -i]));
            }

            H[n] = 20*Math.log10(numerator.divide(denominator).abs());

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

    public double[] normalize(short[] r){
        double[] result = new double[r.length];

        /*for(int i= 0; i<r.length; i++){
            result[i] = r[i]/32768.0;
        }*/

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
        if(factor < 200)
            factor = 1;
        //Log.d("factor",String.valueOf(factor));

        for(int i=0; i<r.length; ++i){
            result[i] = r[i]/factor;
        }

        return result;
    }

    public double[] normalize(double[] r){
        double[] result = new double[r.length];

        /*for(int i= 0; i<r.length; i++){
            result[i] = r[i]/32768.0;
        }*/


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

        for(int i=0; i<r.length; ++i){
            result[i] = r[i]/factor;
        }
        return result;
    }

    public double[] preEmphasis(short[] r){
        double[] result = new double[r.length];
        result[0] = r[0];
        for(int i=1; i< r.length;i++){
            result[i] = r[i] - 0.98*r[i-1];
        }
        return result;
    }

    public double[] hamming(short[] r){
        int N = r.length;
        double[] result = new double[N];
        for(int i=1; i<N -1; ++i){
            double h = 0.54 - 0.46 * Math.cos(2* Math.PI * (double)i / (N -1));
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
