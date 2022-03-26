#include <smpphat/signal.h>

mics_obj * mics_construct(const char * micarray_str) {

    mics_obj * obj;
    unsigned int channel_index;

    // Allocate memory for object
    obj = (mics_obj *) malloc(sizeof(mics_obj));

    // Load mic array configuration
    if (strcmp(micarray_str, "respeaker_usb") == 0) {
        obj->channels_count = MICARRAY_RESPEAKER_USB_COUNT;
        obj->xyzs = (float *) malloc(sizeof(float) * MICARRAY_RESPEAKER_USB_COUNT * 3);
        memcpy(obj->xyzs, MICARRAY_RESPEAKER_USB_XYZ, sizeof(float) * MICARRAY_RESPEAKER_USB_COUNT * 3);
    }
    if (strcmp(micarray_str, "respeaker_core") == 0) {
        obj->channels_count = MICARRAY_RESPEAKER_CORE_COUNT;
        obj->xyzs = (float *) malloc(sizeof(float) * MICARRAY_RESPEAKER_CORE_COUNT * 3);
        memcpy(obj->xyzs, MICARRAY_RESPEAKER_CORE_XYZ, sizeof(float) * MICARRAY_RESPEAKER_CORE_COUNT * 3);
    }
    if (strcmp(micarray_str, "minidsp_uma") == 0) {
        obj->channels_count = MICARRAY_MINIDSP_UMA_COUNT;
        obj->xyzs = (float *) malloc(sizeof(float) * MICARRAY_MINIDSP_UMA_COUNT * 3);
        memcpy(obj->xyzs, MICARRAY_MINIDSP_UMA_XYZ, sizeof(float) * MICARRAY_MINIDSP_UMA_COUNT * 3);
    }
    if (strcmp(micarray_str, "matrix_creator") == 0) {
        obj->channels_count = MICARRAY_MATRIX_CREATOR_COUNT;
        obj->xyzs = (float *) malloc(sizeof(float) * MICARRAY_MATRIX_CREATOR_COUNT * 3);
        memcpy(obj->xyzs, MICARRAY_MATRIX_CREATOR_XYZ, sizeof(float) * MICARRAY_MATRIX_CREATOR_COUNT * 3);
    }

    // Return pointer to object
    return obj;

}

void mics_destroy(mics_obj * obj) {

    // Deallocate memory for xyz positions
    free((void *) obj->xyzs);

    // Deallocate memory for object
    free((void *) obj);

}

void mics_printf(const mics_obj * obj) {

    unsigned int mic_index;

    for (mic_index = 0; mic_index < obj->channels_count; mic_index++) {

        printf("[%02u]: (%+1.6f, %+1.6f, %+1.6f)\n", mic_index, obj->xyzs[mic_index*3+0], obj->xyzs[mic_index*3+1], obj->xyzs[mic_index*3+2]);

    }

}

pdoas_obj * pdoas_construct(const unsigned int levels_count, const char * shape_str) {

    pdoas_obj * obj;

    // Allocate memory for object
    obj = (pdoas_obj *) malloc(sizeof(pdoas_obj));

    // If full sphere, load the corresponding points
    if (strcmp(shape_str, "full") == 0) {
        
        switch(levels_count) {
        
            case 0:
                obj->points_count = SPHR_FULL_0_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_FULL_0_COUNT * 3);
                memcpy(obj->xyzs, SPHR_FULL_0_DATA, sizeof(float) * SPHR_FULL_0_COUNT * 3);
                break;
        
            case 1:
                obj->points_count = SPHR_FULL_1_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_FULL_1_COUNT * 3);
                memcpy(obj->xyzs, SPHR_FULL_1_DATA, sizeof(float) * SPHR_FULL_1_COUNT * 3);
                break;
        
            case 2:
                obj->points_count = SPHR_FULL_2_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_FULL_2_COUNT * 3);
                memcpy(obj->xyzs, SPHR_FULL_2_DATA, sizeof(float) * SPHR_FULL_2_COUNT * 3);
                break;
        
            case 3:
                obj->points_count = SPHR_FULL_3_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_FULL_3_COUNT * 3);
                memcpy(obj->xyzs, SPHR_FULL_3_DATA, sizeof(float) * SPHR_FULL_3_COUNT * 3);
                break;
        
            case 4:
                obj->points_count = SPHR_FULL_4_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_FULL_4_COUNT * 3);
                memcpy(obj->xyzs, SPHR_FULL_4_DATA, sizeof(float) * SPHR_FULL_4_COUNT * 3);
                break;

        }
    }
    // If half sphere, load the corresponding points
    if (strcmp(shape_str, "half") == 0) {
        
        switch(levels_count) {
        
            case 0:
                obj->points_count = SPHR_HALF_0_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_HALF_0_COUNT * 3);
                memcpy(obj->xyzs, SPHR_HALF_0_DATA, sizeof(float) * SPHR_HALF_0_COUNT * 3);
                break;
        
            case 1:
                obj->points_count = SPHR_HALF_1_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_HALF_1_COUNT * 3);
                memcpy(obj->xyzs, SPHR_HALF_1_DATA, sizeof(float) * SPHR_HALF_1_COUNT * 3);
                break;
        
            case 2:
                obj->points_count = SPHR_HALF_2_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_HALF_2_COUNT * 3);
                memcpy(obj->xyzs, SPHR_HALF_2_DATA, sizeof(float) * SPHR_HALF_2_COUNT * 3);
                break;
        
            case 3:
                obj->points_count = SPHR_HALF_3_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_HALF_3_COUNT * 3);
                memcpy(obj->xyzs, SPHR_HALF_3_DATA, sizeof(float) * SPHR_HALF_3_COUNT * 3);
                break;
        
            case 4:
                obj->points_count = SPHR_HALF_4_COUNT;
                obj->xyzs = (float *) malloc(sizeof(float) * SPHR_HALF_4_COUNT * 3);
                memcpy(obj->xyzs, SPHR_HALF_4_DATA, sizeof(float) * SPHR_HALF_4_COUNT * 3);
                break;

        }
    }

    // Return pointer to object
    return obj;

}

void pdoas_destroy(pdoas_obj * obj) {

    // Deallocate memory for xyz positions
    free((void *) obj->xyzs);

    // Deallocate memory for object
    free((void *) obj);

}

void pdoas_printf(const pdoas_obj * obj) {

    unsigned int point_index;

    for (point_index = 0; point_index < obj->points_count; point_index++) {

        printf("[%02u]: (%+1.6f, %+1.6f, %+1.6f)\n", point_index, obj->xyzs[point_index*3+0], obj->xyzs[point_index*3+1], obj->xyzs[point_index*3+2]);

    }

}

ldoas_obj * ldoas_construct(const unsigned int points_count) {

    ldoas_obj * obj;

    // Allocate memory for object
    obj = (ldoas_obj *) malloc(sizeof(ldoas_obj));

    // Create room to store doas
    obj->points_count = points_count;
    obj->xyzs = (float *) malloc(sizeof(float) * points_count * 3);
    obj->es = (float *) malloc(sizeof(float) * points_count);

    // Return pointer to object
    return obj;

}

void ldoas_destroy(ldoas_obj * obj) {

    // Deallocate memory for doas
    free((void *) obj->xyzs);

    // Deallocate memory for energies
    free((void *) obj->es);

    // Deallocate memory for object
    free((void *) obj);

}

hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size) {

    hops_obj * obj;
    unsigned int channel_index;

    // Allocate memory for object
    obj = (hops_obj *) malloc(sizeof(hops_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->hop_size = hop_size;

    // Create an array, with a frame with hop_size samples for each channel
    obj->samples = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->samples[channel_index] = (float *) malloc(sizeof(float) * hop_size);
        memset(obj->samples[channel_index], 0x00, sizeof(float) * hop_size);
    }

    // Return pointer to object
    return obj;

}

void hops_destroy(hops_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for frames
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->samples[channel_index]);
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

void hops_printf(const hops_obj * obj) {

    unsigned int channel_index;
    unsigned int sample_index;

    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        
        for (sample_index = 0; sample_index < obj->hop_size; sample_index++){

            printf("[%02u]-(%03u): %1.3f\n", channel_index, sample_index, obj->samples[channel_index][sample_index]);

        }

    }

}

freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size) {

    freqs_obj * obj;
    unsigned int channel_index;

    // Allocate memory for object
    obj = (freqs_obj *) malloc(sizeof(freqs_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;

    // Create an array, with a frame of frame_size/2+1 complex numbers for each channel
    obj->samples = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->samples[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
        memset(obj->samples[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
    }

    // Return pointer to object
    return obj;

}

void freqs_destroy(freqs_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for frames
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->samples[channel_index]);
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

void freqs_printf(const freqs_obj * obj) {

    unsigned int channel_index;
    unsigned int bin_index;

    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        
        for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

            printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", channel_index, bin_index, obj->samples[channel_index][bin_index*2+0], obj->samples[channel_index][bin_index*2+1]);

        }

    }

}

covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size) {

    covs_obj * obj;
    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    // Allocate memory for object
    obj = (covs_obj *) malloc(sizeof(covs_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;

    // For each pair of microphones, create a frame with (frame_size/2+1) complex numbers
    pair_index = 0;
    obj->samples = (float **) malloc(sizeof(float *) * channels_count * (channels_count-1) / 2);
    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
            obj->samples[pair_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
            memset(obj->samples[pair_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
            pair_index++;
        }
    }

    // Return pointer to object
    return obj;

}

void covs_destroy(covs_obj * obj) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    // Deallocate memory for frames
    pair_index = 0;
    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
            free((void *) obj->samples[pair_index]);
            pair_index++;
        }
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

void covs_printf(const covs_obj * obj) {

    unsigned int pair_index;
    unsigned int bin_index;

    for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count-1)/2); pair_index++) {
        
        for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

            printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", pair_index, bin_index, obj->samples[pair_index][bin_index*2+0], obj->samples[pair_index][bin_index*2+1]);

        }

    }

}