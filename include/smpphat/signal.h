#ifndef __SMPPHAT_SIGNAL
#define __SMPPHAT_SIGNAL

    #include <stdlib.h>
    #include <string.h>
    #include <stdio.h>

    #include "const.h"

    //
    // Signal that contains the positions of the microphones
    //
    typedef struct mics_obj {

        // Number of channels
        unsigned int channels_count;
        // Pointer to an array of xyz coordinates
        float * xyzs;

    } mics_obj;

    //
    // Signal that holds all the potential directions of arrivals
    //
    typedef struct pdoas_obj {

        // Number of points
        unsigned int points_count;
        // Pointer to an array of xyz coordinates
        float * xyzs;

    } pdoas_obj;

    //
    // Signal that holds the localized directions of arrivals
    //
    typedef struct ldoas_obj {

        // Number of points
        unsigned int points_count;
        // Pointer to an array of xyz coordinates
        float * xyzs;
        // Pointer to corresponding energy levels
        float * es;

    } ldoas_obj;

    //
    // Signal that contains chunk of samples to fill frames
    //
    typedef struct hops_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int hop_size;
        // Pointer to an array of arrays of samples
        float ** samples;

    } hops_obj;

    //
    // Signal that contains the frames in the frequency domain
    //
    typedef struct freqs_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int frame_size;
        // Pointer to an array of arrays (size: channels_count) of frequency bins (size: frame_size / 2 + 1)
        float ** samples;

    } freqs_obj;

    //
    // Signal that contains the frames of the spatial covariance matrices
    //
    typedef struct covs_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int frame_size;
        // Pointer to an array of arrays (size: channel_count * (channel_count-1) / 2) of cross-spectra (size: frame_size / 2 + 1)
        float ** samples;

    } covs_obj;

    //
    // Construct the mics object
    //
    // micarray_str             Name of the microphone array
    //
    // (return)                 Pointer to the mics object
    //
    mics_obj * mics_construct(const char * micarray_str);

    //
    // Destroy the mics object
    //
    // obj                      Pointer to the mics object
    //
    // (return)                 None
    //
    void mics_destroy(mics_obj * obj);

    //
    // Print the content of the mics object
    //
    // obj                      Pointer to the mics object
    //
    // (return)                 None, but prints content in console
    //
    void mics_printf(const mics_obj * obj);

    //
    // Construct the pdoas object
    //
    // levels_count             Number of refine levels for icosahedron (between 0 and 4)
    // shape_str                Either 'full' for full sphere, or 'half' for half sphere
    //
    // (return)                 Pointer to the pdoas object
    //
    pdoas_obj * pdoas_construct(const unsigned int levels_count, const char * shape_str);

    //
    // Destroy the pdoas object
    //
    // obj                      Pointer to the pdoas object
    //
    // (return)                 None
    //
    void pdoas_destroy(pdoas_obj * obj);

    //
    // Print the content of the pdoas object
    //
    // obj                      Pointer to the pdoas object
    //
    // (return)                 None, but prints content in console
    //
    void pdoas_printf(const pdoas_obj * obj);

    //
    // Construct the ldoas object
    //
    // points_count             Number of directions of arrival
    // 
    // (return)                 Pointer to the ldoas object
    //
    ldoas_obj * ldoas_construct(const unsigned int points_count);

    //
    // Destroy the ldoas object
    // 
    // obj                      Pointer to the ldoas object
    //
    // (return)                 None
    //
    void ldoas_destroy(ldoas_obj * obj);

    //
    // Construct the hops object
    //
    // channels_count           Number of channels
    // hop_size                 Number of samples per channel
    //
    // (return)                 Pointer to the hops object
    //
    hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size);

    //
    // Destroy the hops object
    //
    // obj                      Pointer to the hops object
    //
    // (return)                 None
    //
    void hops_destroy(hops_obj * obj);

    //
    // Print the content of the hops object
    //
    // obj                      Pointer to the hops object
    //
    // (return)                 None, but prints content in console
    //
    void hops_printf(const hops_obj * obj);

    //
    // Construct the freqs object
    //
    // channels_count           Number of channels
    // frame_size               Number of samples per channel in each frame in the time-domain
    //
    // (return)                 Pointer to the freqs object
    //
    freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size);

    //
    // Destroy the freqs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None
    //
    void freqs_destroy(freqs_obj * obj);

    //
    // Print the content of the freqs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None, but prints content in console
    //
    void freqs_printf(const freqs_obj * obj);

    //
    // Construct the covs object
    //
    // channels_count           Number of channels
    // frame_size               Number of samples per channel in each frame in the time-domain
    //
    // (return)                 Pointer to the covs object
    //
    covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size);

    //
    // Destroy the covs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None
    //
    void covs_destroy(covs_obj * obj);

    //
    // Print the content of the covs object
    //
    // obj                      Pointer to the covs object
    //
    // (return)                 None, but prints content in console
    //
    void covs_printf(const covs_obj * obj);


#endif