package com.example.ryosuke.aasmp;

import android.Manifest;
import android.content.pm.PackageManager;
import android.os.Bundle;
import android.support.annotation.NonNull;
import android.support.design.widget.Snackbar;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;

import com.newventuresoftware.waveform.WaveformView;

import java.util.Arrays;

public class MainActivity extends AppCompatActivity {

    private final int LPC_ORDER = 12;

    private int REQUEST_PERMISSION = 100;
    private RecordingThread mRecordingThread;

    private WaveformView mRealtimeWaveformView;
    private Lpc mLpc;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        mRealtimeWaveformView = (WaveformView) findViewById(R.id.waveformView);
        mRecordingThread = new RecordingThread(new AudioDataReceivedListener() {
            @Override
            public void onAudioDataReceived(short[] data) {

                int center = data.length/2;
                double cuttime = 0.04;
                int SAMPLE_RATE = RecordingThread.SAMPLE_RATE;
                short[] s = Arrays.copyOfRange(data, (int)(center- cuttime*SAMPLE_RATE/2), (int)(center + cuttime*SAMPLE_RATE/2));
                double df = (double)SAMPLE_RATE/s.length;

                double[] hamming_result = mLpc.normalize(mLpc.hamming(s));
                double[] lpc_result = mLpc.normalize(mLpc.lpc(hamming_result, LPC_ORDER, df));
                short[] reNewdata = mLpc.toShort(lpc_result);
                short[] halfData = Arrays.copyOfRange(reNewdata, 0, reNewdata.length/2);

                Lpc.Formant formant_result = mLpc.formant(lpc_result, df);
                double f1 = formant_result.first;
                double f2 = formant_result.second;

                Log.d("Formant", "f1 is" + f1 + ". f2 is"  + f2);



                mRealtimeWaveformView.setSamples(halfData);
            }
        });

        Button button = (Button) findViewById(R.id.button);
        button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if(!mRecordingThread.recording()){
                    startAudioRecordingSafe();
                } else {
                    mRecordingThread.stopRecording();
                }
            }
        });

        mLpc = new Lpc();



    }

    @Override
    protected void onStop() {
        super.onStop();

        mRecordingThread.stopRecording();
    }

    private void startAudioRecordingSafe() {
        if (ContextCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO)
                == PackageManager.PERMISSION_GRANTED) {
            mRecordingThread.startRecording();
        } else {
            requestMicrophonePermission();
        }
    }

    private void requestMicrophonePermission(){
        if(ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.RECORD_AUDIO) && ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.RECORD_AUDIO) ){
            Snackbar.make(mRealtimeWaveformView,"Microphone access is required in order to record audio",
                    Snackbar.LENGTH_INDEFINITE).setAction("OK", new View.OnClickListener(){
                @Override
                public void onClick(View v){
                    ActivityCompat.requestPermissions(MainActivity.this, new String[]{
                            Manifest.permission.RECORD_AUDIO,Manifest.permission.MODIFY_AUDIO_SETTINGS}, REQUEST_PERMISSION);
                    }

            }).show();

        } else {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{
                    Manifest.permission.RECORD_AUDIO, Manifest.permission.MODIFY_AUDIO_SETTINGS},
                    REQUEST_PERMISSION);
        }
    }

    @Override
    public void onRequestPermissionsResult(int requestCode, @NonNull String[] permissions,
                                         @NonNull int[] grantResults){
        if(requestCode == REQUEST_PERMISSION && grantResults.length > 0 &&
                grantResults[0] == PackageManager.PERMISSION_GRANTED){

        }

    }
}
