#ifndef __SMPPHAT_SYSTEM
#define __SMPPHAT_SYSTEM

    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <math.h>
    #include <fftw3.h>
    #include "signal.h"
    #include "const.h"

    //
    // Wave file header
    //
    // This structure holds the first 44 bytes used as a header
    // for the WAVE format.
    //
    typedef struct wav_header {

        // Chunk descriptor (corresponds to characters "RIFF")
        char chunk_id[4];
        // Size of the file in bytes - 8
        unsigned int chunk_size;
        // Contains the characters "WAVE"
        char format[4];
        // Contains the characters "fmt "
        char subchunk1_id[4];
        // Corresponds to 16 for PCM
        unsigned int subchunk1_size;
        // PCM = 1 for linear quantization
        unsigned short audio_format;
        // Number of channels
        unsigned short num_channels;
        // Sample rate in samples/sec
        unsigned int sample_rate;
        // Equals sample_rate * num_channels * bits_per_sample / 8
        unsigned int byte_rate;
        // Equals num_channels * bits_per_sample / 8
        unsigned short block_align;
        // Number of bits per sample
        unsigned short bits_per_sample;
        // Contains the characters "data"
        char subchunk2_id[4];
        // Equals num_samples * num_channels * bits_per_sample / 8
        unsigned int subchunk2_size;

    } wav_header;

    //
    // Wave object that handles reading a WAVE file
    //
    typedef struct wav_obj {

        // Sample rate in samples/sec
        unsigned int sample_rate;
        // Number of channels
        unsigned int num_channels;
        // Number of bits per sample
        unsigned int bits_per_sample;
        // Number of samples read per channel each time
        unsigned int hop_size;

        // File pointer to read from the wave file
        FILE * file_pointer;
        // Buffer to load the samples (here we assume 16 bits per sample)
        short * buffer;

    } wav_obj;

    //
    // STFT object to perform Short-Time Fourier Transform
    //
    typedef struct stft_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per frame to perform FFT
        unsigned int frame_size;
        // Hop size between frames in samples
        unsigned int hop_size;

        // This pointer refers to an array with the window's coefficients
        // (size = frame_size)
        float * window;
        // This pointer refers to an array of frames, each one containing an array of samples
        // (size = num_channels, and for each frame, size = frame_size)
        float ** frames;

        // This is the FFTW plan
        fftwf_plan fft;
        // This array points to the time-domain samples to perform the FFT
        // (size = frame_size)
        float * frame_real;
        // This array points to the frequency-domain samples obtained after FFT
        // (size = frame_size/2+1)
        fftwf_complex * frame_complex;

    } stft_obj;    

    //
    // Spatial Covariance Matrix object to estimate the cross-spectra and
    // perform phase transform to normalize the amplitude in the frequency domain
    //
    typedef struct scmphat_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples in the time-domain prior to the STFT        
        unsigned int frame_size;
        // Smoothing parameter to average over time
        // (alpha is between [0,1], and alpha=1 means no smoothing at all)        
        float alpha;
        // Method to compute the result: if set to 'c', it will normalize each channel
        // before computing the cross-spectrum, and if set to 'p' it will compute the
        // cross-spectrum and then normalize
        char method;

        // Array of frames in the frequency domain with the cross-spectra
        float ** cross_spectrum;

    } scmphat_obj;

    //
    // Steered Response Power with Merged Pairs
    //
    typedef struct smp_obj {

        // Number of samples in the time-domain prior to the STFT
        unsigned int frame_size;
        // Interpolation rate for the IFFT (r=1, r=2, r=4, ...)
        unsigned int interpolation_rate;
        
        // Speed of sound (m/sec)
        float sound_speed;
        // Sample rate
        unsigned int sample_rate;
        // Epsilon value
        float epsilon;

        // Number of potential DoAs
        unsigned int pdoas_count;
        // Coordinates of potential DoAs
        float * pdoas_xyzs;
        // Number of pairs of microphones
        unsigned int pairs_count;
        // Vectors of the difference between mics positions for each pair
        float * pairs_xyzs;

        // Number of groups
        unsigned int groups_count;
        // Array of group indexes
        unsigned int * groups;
        // Array of polarities
        float * polarities;
        // Array of TDoAs
        float * tdoas;
        // Range of TDoAs
        unsigned int * ranges;
        // Lookup for TDoAs
        unsigned int * lookups;

        // These are the FFTW plans
        fftwf_plan * iffts;
        // These arrays point to the frequency-domain samples obtained after FFT
        // (size = frame_size*interpolation_rate/2+1)
        fftwf_complex ** frames_complex;
        // These arrays point to the time-domain samples to perform the FFT
        // (size = frame_size)
        float ** frames_real;

        // Cross-correlation results
        float * cross_correlation;

    } smp_obj;

    //
    // Steered Response Power
    //
    typedef struct srp_obj {

        // Number of samples in the time-domain prior to the STFT
        unsigned int frame_size;
        // Interpolation rate for the IFFT (r=1, r=2, r=4, ...)
        unsigned int interpolation_rate;
        
        // Speed of sound (m/sec)
        float sound_speed;
        // Sample rate
        unsigned int sample_rate;

        // Number of potential DoAs
        unsigned int pdoas_count;
        // Coordinates of potential DoAs
        float * pdoas_xyzs;
        // Number of pairs of microphones
        unsigned int pairs_count;
        // Vectors of the difference between mics positions for each pair
        float * pairs_xyzs;

        // Array of TDoAs
        float * tdoas;
        // Range of TDoAs
        unsigned int * ranges;
        // Lookup for TDoAs
        unsigned int * lookups;

        // This is the FFTW plan
        fftwf_plan ifft;
        // This array points to the frequency-domain samples obtained after FFT
        // (size = frame_size*interpolation_rate/2+1)
        fftwf_complex * frame_complex;
        // This array points to the time-domain samples to perform the FFT
        // (size = frame_size)
        float * frame_real;

        // Cross-correlation results
        float * cross_correlation;        

    } srp_obj;

    //
    // CSV object to write results to file
    //
    typedef struct csv_obj {

        // File pointer
        FILE * file_pointer;

        // Number of channels
        unsigned int channels_count;
        // The frame index
        unsigned int frame_index;

    } csv_obj;

    //
    // Construct the wave object
    //
    // file_name                Path of the wave file
    // hop_size                 Number of samples per channel to read at once
    //
    // (return)                 Pointer to the wave object
    //
    wav_obj * wav_construct(const char * file_name, const unsigned int hop_size);

    //
    // Destroy the wave object
    //
    // obj                      Pointer to the wave object
    //
    // (return)                 None
    //
    void wav_destroy(wav_obj * obj);

    //
    // Read samples from the wave object
    //
    // obj                  Pointer to the wave object
    // hops                 Pointer to the hops object that receives samples
    //
    // (return)             Returns -1 if the end of file reached, 0 otherwise
    //
    int wav_read(wav_obj * obj, hops_obj * hops);

    //
    // Construct the stft object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // hop_size             Number of samples between frame
    //
    // (return)             Pointer of the stft object
    //
    stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size);

    //
    // Destroy the stft object
    //
    // obj                  Pointer to the stft object
    //
    // (return)             None
    //
    void stft_destroy(stft_obj * obj);

    //
    // Compute the Short-Time Fourier Transform
    //
    // obj                  Pointer to the stft object
    // hops                 Pointer to the hops object with new samples
    // freqs                Pointer to the freqs object with the results
    //
    // (return)             Returns 0 if no error
    //
    int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs);

    // 
    // Construct the scmphat object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // alpha                Smoothing factor
    //
    // (return)             Pointer to the scmphat object
    //
    scmphat_obj * scmphat_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha);

    //
    // Destroy the scmphat object
    //
    // obj                  Pointer to the scmphat object
    //
    // (return)             None
    //
    void scmphat_destroy(scmphat_obj * obj);

    //
    // Compute the Spatial Covariance Matrix and normalize
    //
    // obj                  Pointer to the scmphat object
    // freqs                Pointer to the freqs object with stft samples
    // covs                 Pointer to the covs object with scm results
    //
    // (return)             Returns 0 if no error
    //
    int scmphat_call(scmphat_obj * obj, const freqs_obj * freqs, covs_obj * covs);

    //
    // Construct the smp object
    //
    // frame_size           Number of samples per frame
    // interpolation_rate   Interpolation rate for ifft (r=1, r=2, r=4, ...)
    // sound_speed          Speed of sound (in m/sec)
    // sample_rate          Sample rate (in samples/sec)
    // mics                 Pointer to the mics object
    // pdoas                Pointer to the pdoas object (doas to scan)
    // 
    // (return)             Pointer to the smp object
    //
    smp_obj * smp_construct(const unsigned int frame_size, const unsigned int interpolation_rate, const float sound_speed, const unsigned int sample_rate, const mics_obj * mics, const pdoas_obj * pdoas);

    //
    // Destroy the smp object
    //
    // obj                  Pointer to the smp object
    //
    // (return)             None
    //
    void smp_destroy(smp_obj * obj);

    //
    // Scan for directions of arrival with SMP
    //
    // obj                  Pointer to the smp object
    // covs                 Pointer to the covs object
    // ldoas                Pointer to the ldoas object with ldoas results
    //
    // (return)             Returns 0 if no error
    //
    int smp_call(smp_obj * obj, const covs_obj * covs, ldoas_obj * ldoas);

    //
    // Construct the srp object
    //
    // frame_size           Number of samples per frame
    // interpolation_rate   Interpolation rate for ifft (r=1, r=2, r=4, ...)
    // sound_speed          Speed of sound (in m/sec)
    // sample_rate          Sample rate (in samples/sec)
    // mics                 Pointer to the mics object
    // pdoas                Pointer to the pdoas object (doas to scan)
    // 
    // (return)             Pointer to the srp object
    //
    srp_obj * srp_construct(const unsigned int frame_size, const unsigned int interpolation_rate, const float sound_speed, const unsigned int sample_rate, const mics_obj * mics, const pdoas_obj * pdoas);

    //
    // Destroy the smp object
    //
    // obj                  Pointer to the srp object
    //
    // (return)             None
    //
    void srp_destroy(srp_obj * obj);

    //
    // Scan for directions of arrival with SRP
    //
    // obj                  Pointer to the smp object
    // covs                 Pointer to the covs object
    // ldoas                Pointer to the ldoas object with ldoas results
    //
    // (return)             Returns 0 if no error
    //
    int srp_call(srp_obj * obj, const covs_obj * covs, ldoas_obj * ldoas);

    //
    // Construct the csv object
    //
    // file_name            Path of the file name to write the CSV content
    // channels_count       Number of channels
    //
    // (return)             Pointer to the csv object
    // 
    csv_obj * csv_construct(const char * file_name, const unsigned int channels_count);

    //
    // Destroy the csv object
    //
    // obj                  Pointer to the csv object
    //
    // (return)             None
    //
    void csv_destroy(csv_obj * obj);

    //
    // Write the TDoAs to the CSV file
    //
    // obj                  Pointer to the csv object
    // taus                 Pointer to the ldoas object which contains the DoAs
    //
    // (return)             Returns 0 if no error
    //
    int csv_write(csv_obj * obj, ldoas_obj * ldoas);


#endif