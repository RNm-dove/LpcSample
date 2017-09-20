package com.example.ryosuke.aasmp;

import android.util.Log;

import org.apache.commons.math3.complex.*;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;

/**
 * Created by ryosuke on 17/08/26.
 */

public class Lpc {

    public Lpc() {

    }

    public double[] lpc(double[] input, int order, int fs){

        int N = input.length;

        //自己相関関数
        double[] r = new double[N];

        for(int l=0; l <= order; l++ ){
            r[l] = 0.0;
            for(int n=0; n < N-l; n++){
                r[l] += input[n] * input[n + l];
            }

        }

        //levinson-durbin
        double a[] = new double[order + 1];   //a[0]は１で固定のためlpcOrder個の係数を得るためには+1が必要
        double e[] = new double[order + 1];
        Arrays.fill(a, 0.0);
        Arrays.fill(e, 0.0);

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

        Complex_F64[] rts = findRoots(a);
        ArrayList<Double> arglist = new ArrayList<Double>();
        ArrayList<Complex_F64> rtslist = new ArrayList<Complex_F64>();
        for(int i=0; i< rts.length; i++){

            if(rts[i].getImaginary()>0.0){ //虚部が負の場合を除く
                double arg = Math.atan2(rts[i].getImaginary(), rts[i].getReal());
                arglist.add(arg);
                rtslist.add(rts[i]);
            }
        }
        Integer[] sort_index = argsort(arglist);

        //iterator を用いて角度の大きさの順番を取得。昇順にし、またそれを用いて元の rtslistの順番も並び替える。

        Iterator<Integer> itr = Arrays.asList(sort_index).iterator();

        ArrayList<Double> freqs = new ArrayList<Double>();
        ArrayList<Complex_F64> bw = new ArrayList<Complex_F64>();


        while(itr.hasNext()){
            int num = itr.next();
            freqs.add(arglist.get(num)); //角度昇順
            bw.add(rtslist.get(num));  //rootを昇順
        }

        double[] bwlist = new double[rtslist.size()];

        //角度からフォルマントを計算。　またフォルマント帯域も計算・

        for(int i = 0; i<arglist.size(); i++){
            double arg1 = freqs.get(i)*(fs / (2* Math.PI));
            freqs.set(i, arg1) ;
            double arg2 =  -1.0 /2 * (fs / (2 *Math.PI)) * Math.log(bw.get(i).getMagnitude());
            bwlist[i] = arg2;
        }

        double[] formant = new double[freqs.size()];

        for(int i=0; i<freqs.size(); i++){
            if(freqs.get(i) > 90 && bwlist[i] < 400){
                formant[i] = freqs.get(i);
            }
        }



        return formant;



    }

    public static Integer[] argsort(final ArrayList<Double> a) {
        return argsort(a, true);
    }

    public static Integer[] argsort(final ArrayList<Double> a, final boolean ascending) {
        Integer[] indexes = new Integer[a.size()];
        for (int i = 0; i < indexes.length; i++) {
            indexes[i] = i;
        }
        Arrays.sort(indexes, new Comparator<Integer>() {
            @Override
            public int compare(final Integer i1, final Integer i2) {
                return (ascending ? 1 : -1) * Double.compare(a.get(i1), a.get(i2));
            }
        });
        /*int[] int_value = new int[indexes.length];
        for(int i = 0; i < indexes.length; i++){
            int_value[i] = indexes[i].intValue();
        }*/
        return indexes;
    }


    public static Complex_F64[] findRoots(double... coefficients){
        int N = coefficients.length -1;

        //Construct the companion matrix
        DMatrixRMaj c =new DMatrixRMaj(N, N);

        double a = coefficients[N];
        for ( int i = 0; i < N; i++){
            c.set(i, N-1, -coefficients[i]/a);
        }
        for( int i= 1; i < N; i++){
            c.set(i,i-1,1);
        }

        // use generalized eigenvalue decomposition to find the roots
        EigenDecomposition_F64<DMatrixRMaj> evd = DecompositionFactory_DDRM.eig(N, false);

        evd.decompose(c);

        Complex_F64[] roots = new Complex_F64[N];

        for( int i = 0; i<N ; i++){
            roots[i] = evd.getEigenvalue(i);

        }
        return roots;
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

            H[n] = 20.0*Math.log10(numerator.divide(denominator).abs());

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

    public double[] hamming(double[] r){
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
        return Math.sqrt(v);
    }

    public double volume(short[] r){
        double v = 0;
        for(double x: r){
            v += x*x;
        }
        v /= r.length;
        return Math.sqrt(v);
    }

    public boolean isAudible(short[] data ){
        double rms = volume(data);
        return (rms > 198 && 5600 >rms);
    }

    public String vowel(double f1, double f2)
    {
        if (f1 > 600 && f1 < 1400 && f2 > 900  && f2 < 2000) return "あ";
        if (f1 > 100 && f1 < 410  && f2 > 1900 && f2 < 3500) return "い";
        if (f1 > 100 && f1 < 300  && f2 > 1100 && f2 < 2000) return "う";
        if (f1 > 400 && f1 < 800  && f2 > 1700 && f2 < 3000) return "え";
        if (f1 > 300 && f1 < 900  && f2 > 500  && f2 < 1300) return "お";
        return "-";
    }


}
