package com.example.ryosuke.aasmp;

import android.Manifest;
import android.content.pm.PackageManager;
import android.os.Bundle;
import android.os.Handler;
import android.support.annotation.NonNull;
import android.support.design.widget.Snackbar;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.LinearLayout;
import android.widget.TextView;

import java.util.Arrays;

public class MainActivity extends AppCompatActivity {

    private final int LPC_ORDER = 53;

    private int REQUEST_PERMISSION = 100;
    private RecordingThread mRecordingThread;

    private LinearLayout mGroundView;
    private TextView mTextView;

    private Lpc mLpc;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        setContentView(R.layout.activity_main);

        mGroundView = (LinearLayout) findViewById(R.id.ground_view);
        mTextView = (TextView) findViewById(R.id.text_view);

        mLpc = new Lpc();

        final Handler handler = new Handler();

        mRecordingThread = new RecordingThread(new AudioDataReceivedListener() {

            @Override
            public void onAudioDataReceived(short[] data) {


                if (mLpc == null) {
                    mLpc = new Lpc();
                }

                if (mLpc.isAudible(data)) {
                    int center = data.length / 2;
                    double cuttime = 0.04;
                    int SAMPLE_RATE = RecordingThread.getSampleRate();
                    short[] s = Arrays.copyOfRange(data, (int) (center - cuttime * SAMPLE_RATE / 2), (int) (center + cuttime * SAMPLE_RATE / 2));


                    double df = (double) SAMPLE_RATE / s.length;
                    //double df = SAMPLE_RATE/(double)data.length;

                    double[] double_result = mLpc.toDouble(s);
                    double[] pre_result = mLpc.preEmphasis(double_result);
                    double[] hamming_result = mLpc.hamming(pre_result);
                    double[] formant = mLpc.lpc(hamming_result, LPC_ORDER, SAMPLE_RATE, df);

                    double f1 = formant[0];
                    double f2 = formant[1];
                    final String result_string = mLpc.vowel(f1, f2);

                    handler.post(new Runnable() {
                        @Override
                        public void run() {
                            if(!result_string.equals("")){
                                mTextView.setText(result_string);
                            }
                        }
                    });

                    Log.d("Formant", "f1:" + f1 + ",f2:"  + f2 + "," + mLpc.vowel(f1,f2));

                }


            }
        });

        Button button = (Button) findViewById(R.id.button);
        button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (!mRecordingThread.recording()) {
                    startAudioRecordingSafe();
                } else {
                    mRecordingThread.stopRecording();
                }
            }
        });

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

    private void requestMicrophonePermission() {
        if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.RECORD_AUDIO) && ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.RECORD_AUDIO)) {
            Snackbar.make(mGroundView, "Microphone access is required in order to record audio",
                    Snackbar.LENGTH_INDEFINITE).setAction("OK", new View.OnClickListener() {
                @Override
                public void onClick(View v) {
                    ActivityCompat.requestPermissions(MainActivity.this, new String[]{
                            Manifest.permission.RECORD_AUDIO, Manifest.permission.MODIFY_AUDIO_SETTINGS}, REQUEST_PERMISSION);
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
                                           @NonNull int[] grantResults) {
        if (requestCode == REQUEST_PERMISSION && grantResults.length > 0 &&
                grantResults[0] == PackageManager.PERMISSION_GRANTED) {

        }

    }

}


