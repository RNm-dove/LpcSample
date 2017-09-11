package com.example.ryosuke.aasmp;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.os.Process;
import android.util.Log;

/**
 * Created by ryosuke on 17/09/07.
 */

public class RecordingThread {

    private static final String LOG_TAG = RecordingThread.class.getSimpleName();
    private static final int SAMPLE_RATE = 44100;

    private boolean isRecording;
    private Thread mThread;
    private AudioDataReceivedListener mListener;

    private double mDf;

    public RecordingThread(AudioDataReceivedListener listener){mListener = listener;}

    public boolean recording(){return mThread != null;}

    public void startRecording(){
        if(mThread != null)
            return;

        isRecording = true;
        mThread = new Thread(new Runnable() {
            @Override
            public void run() {
                record();
            }
        });
        mThread.start();

    }

    public void stopRecording(){
        if(mThread == null)
            return;

        isRecording = false;
        mThread = null;
    }

    private void record(){
        Log.v(LOG_TAG, "Start");
        Process.setThreadPriority(Process.THREAD_PRIORITY_AUDIO);

        //buffer size in bytes
        int bufferSize = AudioRecord.getMinBufferSize(SAMPLE_RATE,
                AudioFormat.CHANNEL_IN_MONO,
                AudioFormat.ENCODING_PCM_16BIT);

        if(bufferSize == AudioRecord.ERROR || bufferSize == AudioRecord.ERROR_BAD_VALUE){
            bufferSize = SAMPLE_RATE * 2;
        }

        short[] audioBuffer = new short[bufferSize / 2];

        mDf = (double)SAMPLE_RATE/audioBuffer.length;

        AudioRecord record = new AudioRecord(MediaRecorder.AudioSource.DEFAULT,
                SAMPLE_RATE,
                AudioFormat.CHANNEL_IN_MONO,
                AudioFormat.ENCODING_PCM_16BIT,
                bufferSize);

        if(record.getState() != AudioRecord.STATE_INITIALIZED){
            Log.e(LOG_TAG, "Audio Record cant initialize!");
            return;
        }
        record.startRecording();

        Log.v(LOG_TAG, "Start Recording");

        long shortsRead = 0;
        while (isRecording){
            int numberOfShort = record.read(audioBuffer, 0, audioBuffer.length);
            shortsRead += numberOfShort;

            mListener.onAudioDataReceived(audioBuffer, mDf);

        }

        record.stop();
        record.release();

        Log.v(LOG_TAG, String.format("Recording stopped. Samples read: %d", shortsRead));

    }

}
