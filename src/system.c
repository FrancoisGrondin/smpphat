#include <smpphat/system.h>

wav_obj * wav_construct(const char * file_name, const unsigned int hop_size) {

    wav_obj * obj;
    wav_header hdr;
    size_t rtn;

    // Allocate memory for object
    obj = (wav_obj *) malloc(sizeof(wav_obj));

    // Open the wave file to be read
    obj->file_pointer = fopen(file_name, "rb");

    // Load the first 44 bytes that are the header of the file
    rtn = fread(&hdr, sizeof(char), sizeof(wav_header), obj->file_pointer);

    // Copy relevant parameters
    obj->sample_rate = hdr.sample_rate;
    obj->num_channels = hdr.num_channels;
    obj->bits_per_sample = hdr.bits_per_sample;
    obj->hop_size = hop_size;

    // Here we only support signed 16-bit samples
    if (obj->bits_per_sample != 16) {
        printf("Unsupported PCM format.\n");
        exit(-1);
    }

    // Create a buffer of shorts with a total of (hop_size * num_channels) samples
    obj->buffer = (short *) malloc(sizeof(short) * obj->hop_size * obj->num_channels);

    // Return the pointer to the wave object
    return obj;

}

void wav_destroy(wav_obj * obj) {

    // First close the wave file
    fclose(obj->file_pointer);
    // Then free the buffer that held all the samples
    free((void *) obj->buffer);

    // Finally deallocate memory for object
    free((void *) obj);

}

int wav_read(wav_obj * obj, hops_obj * hops) {
    
    size_t rtn;
    unsigned int channel_index;
    unsigned int sample_index;

    // Read from the file a total of 2 * hop_size * num_channels bytes
    rtn = fread(obj->buffer, obj->bits_per_sample/8, obj->hop_size * obj->num_channels, obj->file_pointer);

    // Then load each sample in the buffer (samples are interleaved, and are signed 16-bit)
    // Then are divided by 32768 to be normalized between -1 and 1
    for (sample_index = 0; sample_index < obj->hop_size; sample_index++) {
        for (channel_index = 0; channel_index < obj->num_channels; channel_index++) {
            hops->samples[channel_index][sample_index] = ((float) obj->buffer[sample_index * obj->num_channels + channel_index]) / 32768.0;
        }
    }

    // If the number of elements read is consisten, return 0 (meaning there are still elements to read)
    // otherwise return -1, which means we've reached the end of the file
    return (rtn == obj->hop_size * obj->num_channels) ? 0:-1;

}

stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size) {

    stft_obj * obj;
    unsigned int channel_index;
    unsigned int sample_index;

    // Allocate memory for object
    obj = (stft_obj *) malloc(sizeof(stft_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->hop_size = hop_size;

    // Create a Hann window
    obj->window = (float *) malloc(sizeof(float) * frame_size);
    for (sample_index = 0; sample_index < frame_size; sample_index++) {
        obj->window[sample_index] = powf(sinf((float) sample_index * M_PI / ((float) (frame_size - 1))), 2.0);
    }

    // Allocate memory for each frame. There is one frame of frame_size samples per channel.
    obj->frames = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->frames[channel_index] = (float *) malloc(sizeof(float) * frame_size);
        memset((void *) obj->frames[channel_index], 0x00, sizeof(float) * frame_size);
    }

    // Allocate memory to perform the FFT using FFTW
    obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size);
    obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (frame_size/2+1));
    obj->fft = fftwf_plan_dft_r2c_1d(frame_size, obj->frame_real, obj->frame_complex, FFTW_ESTIMATE);

    // Return the pointer to the stft object
    return obj;

}

void stft_destroy(stft_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for window
    free((void *) obj->window);

    // Deallocate memory for each frame, and then clear list of pointers
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->frames[channel_index]);
    }
    free((void *) obj->frames);

    // Deallocate memory for FFTW
    fftwf_free(obj->frame_real);
    fftwf_free(obj->frame_complex);
    fftwf_destroy_plan(obj->fft);
    fftwf_cleanup();

    // Finally, deallocate memory for the object
    free((void *) obj);

}

int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs) {

    unsigned int channel_index;
    unsigned int sample_index;
    unsigned int bin_index;

    // Load the samples for each channel
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

        //
        // Shift samples to the left:
        // 
        // The initial array 
        // 
        // |AAAAAAAAAAAAAA|BBBBBBBBBBBBBBBBBBBBBBBBBBB|
        //  <--hop_size--> <--(frame_size-hop_size)--> 
        //
        // becomes       
        // 
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|00000000000000|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        memmove(&(obj->frames[channel_index][0]), &(obj->frames[channel_index][obj->hop_size]), sizeof(float) * (obj->frame_size - obj->hop_size));

        //
        // Then fill the right side with new samples
        //
        // So the array 
        // 
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|00000000000000|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        // becomes
        //
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|CCCCCCCCCCCCCC|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        memcpy(&(obj->frames[channel_index][obj->frame_size - obj->hop_size]), &(hops->samples[channel_index][0]), sizeof(float) * obj->hop_size);

        // Then loop for each sample in the frame and multiply it by the window, and save result to array to perform FFT
        for (sample_index = 0; sample_index < obj->frame_size; sample_index++) {
            obj->frame_real[sample_index] = obj->window[sample_index] * obj->frames[channel_index][sample_index];
        }

        // Perform FFT (result goes in frame_complex)
        fftwf_execute(obj->fft);

        // Extract the results and copy them to the freqs object
        for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
            freqs->samples[channel_index][bin_index * 2 + 0] = obj->frame_complex[bin_index][0];
            freqs->samples[channel_index][bin_index * 2 + 1] = obj->frame_complex[bin_index][1];
        }

    }

    // By default return 0 when returns without error
    return 0;

}

scmphat_obj * scmphat_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha) {

    scmphat_obj * obj;
    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;
    unsigned int channel_index;

    // Allocate memory for the object
    obj = (scmphat_obj *) malloc(sizeof(scmphat_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->alpha = alpha;

    // If we take only the instantaneous cross-correlation (alpha = 1)
    // and there are more than 2 channels, it is faster to perform PHAT
    // normalization on each channel and then do the pairwise multiplication
    obj->method = 'p';
    if (alpha == 1.0f) {
        if (channels_count > 2) {
            obj->method = 'c';
        }
    }

    // If pairwise approach, then assign the cross spectrum arrays
    if (obj->method == 'p') {

        pair_index = 0;
        obj->cross_spectrum = (float **) malloc(sizeof(float *) * (channels_count * (channels_count-1) / 2));
        for (channel_index1 = 0; channel_index1 < channels_count; channel_index1++) {
            for (channel_index2 = (channel_index1+1); channel_index2 < channels_count; channel_index2++) {
                obj->cross_spectrum[pair_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
                memset((void *) obj->cross_spectrum[pair_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
                pair_index++;
            }
        }

    }
    // If channelwise approach, then assign the cross spectrum arrays
    if (obj->method == 'c') {

        obj->cross_spectrum = (float **) malloc(sizeof(float *) * channels_count);
        for (channel_index = 0; channel_index < channels_count; channel_index++) {
            obj->cross_spectrum[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
            memset((void *) obj->cross_spectrum[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
        }

    }

    // Return pointer to the object
    return obj;

}

void scmphat_destroy(scmphat_obj * obj) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;
    unsigned int channel_index;

    // If pairwise approach, deallocate memory accordingly
    if (obj->method == 'p') {

        pair_index = 0;
        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
            for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
                free((void *) obj->cross_spectrum[pair_index]);
                pair_index++;
            }
        }
        free((void *) obj->cross_spectrum);

    }
    // Or deallocate in the channelwise approach
    if (obj->method == 'c') {

        for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
            free((void *) obj->cross_spectrum[channel_index]);
        }
        free((void *) obj->cross_spectrum);

    }

    // Deallocate memory for object
    free((void *) obj);

}

int scmphat_call(scmphat_obj * obj, const freqs_obj * freqs, covs_obj * covs) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int bin_index;
    unsigned int pair_index;
    unsigned int channel_index;

    float dest_real;
    float dest_imag;
    float dest_magn;
    float src1_real;
    float src1_imag;
    float src2_real;
    float src2_imag;
    float src_real;
    float src_imag;
    float src_magn;

    // If pairwise approach
    if (obj->method == 'p') {

        // Compute the following for each pair of microphones i and j:
        //
        // 1) R_(i,j)(t,f) = (1 - alpha) * R_(i,j)(t-1,f) + alpha * X_i(t,f) * X_j(t,f)^*
        // 2) Rhat_(i,j)(t,f) = R_(i,j)(t,f) / |R_(i,j)(t,f)|

        pair_index = 0;

        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

            for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

                for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                    dest_real = obj->cross_spectrum[pair_index][bin_index*2+0];
                    dest_imag = obj->cross_spectrum[pair_index][bin_index*2+1];

                    src1_real = freqs->samples[channel_index1][bin_index*2+0];
                    src1_imag = freqs->samples[channel_index1][bin_index*2+1];
                    src2_real = freqs->samples[channel_index2][bin_index*2+0];
                    src2_imag = freqs->samples[channel_index2][bin_index*2+1];                

                    dest_real *= (1.0 - obj->alpha);
                    dest_imag *= (1.0 - obj->alpha);

                    dest_real += obj->alpha * (src1_real * src2_real + src1_imag * src2_imag);
                    dest_imag += obj->alpha * (src1_imag * src2_real - src1_real * src2_imag);

                    obj->cross_spectrum[pair_index][bin_index*2+0] = dest_real;
                    obj->cross_spectrum[pair_index][bin_index*2+1] = dest_imag;

                    dest_magn = sqrtf(dest_real * dest_real + dest_imag * dest_imag);

                    dest_real /= (dest_magn + 1e-10);
                    dest_imag /= (dest_magn + 1e-10);

                    covs->samples[pair_index][bin_index*2+0] = dest_real;
                    covs->samples[pair_index][bin_index*2+1] = dest_imag;

                }

                pair_index++;

            }

        }

    }
    // If channelwise approach
    if (obj->method == 'c') {

        // Compute the following
        //
        // Rhat_(i,j)(t,f) = (X_i(t,f)/|X_i(t,f)|) * (X_j(t,f)^*/|X_j(t,f)|)

        for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

            for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                src_real = freqs->samples[channel_index][bin_index*2+0];
                src_imag = freqs->samples[channel_index][bin_index*2+1];
                src_magn = sqrtf(src_real * src_real + src_imag * src_imag);

                src_real /= (src_magn + 1e-10);
                src_imag /= (src_magn + 1e-10);

                obj->cross_spectrum[channel_index][bin_index*2+0] = src_real;
                obj->cross_spectrum[channel_index][bin_index*2+1] = src_imag;

            }

        }

        pair_index = 0;

        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

            for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

                for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                    src1_real = obj->cross_spectrum[channel_index1][bin_index*2+0];
                    src1_imag = obj->cross_spectrum[channel_index1][bin_index*2+1];
                    src2_real = obj->cross_spectrum[channel_index2][bin_index*2+0];
                    src2_imag = obj->cross_spectrum[channel_index2][bin_index*2+1];

                    dest_real = src1_real * src2_real + src1_imag * src2_imag;
                    dest_imag = src1_imag * src2_real - src1_real * src2_imag;

                    covs->samples[pair_index][bin_index*2+0] = dest_real;
                    covs->samples[pair_index][bin_index*2+1] = dest_imag;

                }

                pair_index++;

            }

        }

    }

    // By default return 0 when returns without error
    return 0;

}

smp_obj * smp_construct(const unsigned int frame_size, const unsigned int interpolation_rate, const float sound_speed, const unsigned int sample_rate, const mics_obj * mics, const pdoas_obj * pdoas) {

    smp_obj * obj;

    unsigned int mic_index1;
    unsigned int mic_index2;
    unsigned int pair_index;
    unsigned int pair_index1;
    unsigned int pair_index2;
    float dot_12;
    float diff_1_magn;
    float diff_2_magn;
    unsigned int group_index;
    unsigned int pdoa_index;
    float tdoa;
    unsigned int offset;

    // Allocate memory for the object
    obj = (smp_obj *) malloc(sizeof(smp_obj));

    // Copy relevant parameters
    obj->frame_size = frame_size;
    obj->interpolation_rate = interpolation_rate;
    obj->sound_speed = sound_speed;
    obj->sample_rate = sample_rate;
    obj->epsilon = 1E-5;
    obj->pdoas_count = pdoas->points_count;
    obj->pdoas_xyzs = (float *) malloc(sizeof(float) * pdoas->points_count * 3);
    memcpy(obj->pdoas_xyzs, pdoas->xyzs, sizeof(float) * pdoas->points_count * 3);
    
    // Compute vectors for each pair
    obj->pairs_count = (mics->channels_count * (mics->channels_count-1) / 2);
    obj->pairs_xyzs = (float *) malloc(sizeof(float) * obj->pairs_count * 3);
    pair_index = 0;
    for (mic_index1 = 0; mic_index1 < mics->channels_count; mic_index1++) {
        for (mic_index2 = (mic_index1+1); mic_index2 < mics->channels_count; mic_index2++) {
            obj->pairs_xyzs[pair_index*3+0] = mics->xyzs[mic_index1*3+0] - mics->xyzs[mic_index2*3+0];
            obj->pairs_xyzs[pair_index*3+1] = mics->xyzs[mic_index1*3+1] - mics->xyzs[mic_index2*3+1];
            obj->pairs_xyzs[pair_index*3+2] = mics->xyzs[mic_index1*3+2] - mics->xyzs[mic_index2*3+2];
            pair_index++;
        }
    }

    // Create merging plan
    obj->groups = (unsigned int *) malloc(sizeof(unsigned int) * obj->pairs_count);
    memset(obj->groups, 0x00, sizeof(unsigned int) * obj->pairs_count);
    obj->polarities = (float *) malloc(sizeof(float) * obj->pairs_count);
    memset(obj->polarities, 0x00, sizeof(float) * obj->pairs_count);

    // Merge pairs together
    group_index = 0;
    for (pair_index1 = 0; pair_index1 < obj->pairs_count; pair_index1++) {

        // If this pair is not associated to a group yet
        if (obj->polarities[pair_index1] == 0.0) {

            // Create new group
            obj->groups[pair_index1] = group_index;
            obj->polarities[pair_index1] = +1.0;

            for (pair_index2 = 0; pair_index2 < obj->pairs_count; pair_index2++) {

                // If this pair is not associated to a group yet
                if (obj->polarities[pair_index2] == 0) {

                    // Compute the dot product
                    dot_12 = 0.0;
                    dot_12 += obj->pairs_xyzs[pair_index1*3+0] * obj->pairs_xyzs[pair_index2*3+0];
                    dot_12 += obj->pairs_xyzs[pair_index1*3+1] * obj->pairs_xyzs[pair_index2*3+1];
                    dot_12 += obj->pairs_xyzs[pair_index1*3+2] * obj->pairs_xyzs[pair_index2*3+2];

                    // Compute magnitudes
                    diff_1_magn = 0.0;
                    diff_1_magn += obj->pairs_xyzs[pair_index1*3+0] * obj->pairs_xyzs[pair_index1*3+0];
                    diff_1_magn += obj->pairs_xyzs[pair_index1*3+1] * obj->pairs_xyzs[pair_index1*3+1];
                    diff_1_magn += obj->pairs_xyzs[pair_index1*3+2] * obj->pairs_xyzs[pair_index1*3+2];
                    diff_1_magn = sqrtf(diff_1_magn);
                    diff_2_magn = 0.0;
                    diff_2_magn += obj->pairs_xyzs[pair_index2*3+0] * obj->pairs_xyzs[pair_index2*3+0];
                    diff_2_magn += obj->pairs_xyzs[pair_index2*3+1] * obj->pairs_xyzs[pair_index2*3+1];
                    diff_2_magn += obj->pairs_xyzs[pair_index2*3+2] * obj->pairs_xyzs[pair_index2*3+2];
                    diff_2_magn = sqrtf(diff_2_magn);

                    // Check if parallel
                    if (fabs(fabs(dot_12) - diff_1_magn * diff_2_magn) < obj->epsilon) {

                        // Check if equidistant
                        if (fabs(diff_1_magn - diff_2_magn) < obj->epsilon) {

                            // Associate this pair to current group
                            obj->groups[pair_index2] = group_index;

                            // And save polarity
                            if (dot_12 > 0) {
                                obj->polarities[pair_index2] = +1.0;
                            }
                            else {
                                obj->polarities[pair_index2] = -1.0;
                            }

                        }

                    }


                }

            }

            group_index++;

        }

    }

    // Save the total number of groups for later
    obj->groups_count = group_index;

    // Compute the TDoAs for the reference pair in each group
    obj->tdoas = (float *) malloc(sizeof(float) * obj->groups_count * obj->pdoas_count);

    for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {        

        for (group_index = 0; group_index < obj->groups_count; group_index++) {

            for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {

                // If the current pair is the first occurent in the group
                if (obj->groups[pair_index] == group_index) {

                    // tdoa = (fs/c) * (d.u)
                    tdoa = 0.0;
                    tdoa += obj->pdoas_xyzs[pdoa_index*3+0] * obj->pairs_xyzs[pair_index*3+0];
                    tdoa += obj->pdoas_xyzs[pdoa_index*3+1] * obj->pairs_xyzs[pair_index*3+1];
                    tdoa += obj->pdoas_xyzs[pdoa_index*3+2] * obj->pairs_xyzs[pair_index*3+2];
                    tdoa *= (float) sample_rate;
                    tdoa /= sound_speed;
                    tdoa *= -1.0;

                    obj->tdoas[pdoa_index * obj->groups_count + group_index] = tdoa;

                    break;

                }

            }

        }

    }

    // Compute the range for each pair
    obj->ranges = (unsigned int *) malloc(sizeof(unsigned int) * obj->groups_count); 
    for (group_index = 0; group_index < obj->groups_count; group_index++) {
        tdoa = 0.0;
        for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {
            if (fabs(obj->tdoas[pdoa_index * obj->groups_count + group_index]) > tdoa) {
                tdoa = fabs(obj->tdoas[pdoa_index * obj->groups_count + group_index]);
            }
        }
        obj->ranges[group_index] = ((unsigned int) ceilf(tdoa * (float) interpolation_rate)) * 2 + 1;
    }

    // Compute number of samples in the super cross-correlation vector
    offset = 0;
    for (group_index = 0; group_index < obj->groups_count; group_index++) {
        offset += obj->ranges[group_index];
    }

    // Generate the lookup table
    offset = 0;
    obj->lookups = (unsigned int *) malloc(sizeof(unsigned int) * obj->groups_count * obj->pdoas_count);    
    for (group_index = 0; group_index < obj->groups_count; group_index++) {
        for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) { 
            tdoa = obj->tdoas[pdoa_index * obj->groups_count + group_index];
            obj->lookups[pdoa_index * obj->groups_count + group_index] = offset + (obj->ranges[group_index]-1)/2 + ((unsigned int) roundf(tdoa * (float) interpolation_rate));
        }
        offset += obj->ranges[group_index];
    }

    // Create a cross-correlation supervector
    obj->cross_correlation = (float *) malloc(sizeof(float) * offset);

    // Create the FFT setup
    obj->frames_complex = (fftwf_complex **) malloc(sizeof(fftwf_complex *) * obj->groups_count);
    obj->frames_real = (float **) malloc(sizeof(float *) * obj->groups_count);
    obj->iffts = (fftwf_plan *) malloc(sizeof(fftwf_plan) * obj->groups_count);

    for (group_index = 0; group_index < obj->groups_count; group_index++) {
        obj->frames_complex[group_index] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * ((frame_size * interpolation_rate)/2+1));
        obj->frames_real[group_index] = (float *) fftwf_malloc(sizeof(float) * frame_size * interpolation_rate);
        obj->iffts[group_index] = fftwf_plan_dft_c2r_1d(frame_size * interpolation_rate, obj->frames_complex[group_index], obj->frames_real[group_index], FFTW_ESTIMATE);
    }

    return obj;

}

void smp_destroy(smp_obj * obj) {

    unsigned int group_index;

    // Free memory allocated for doas
    free((void *) obj->pdoas_xyzs);

    // Free memory allocated to store pairs differences
    free((void *) obj->pairs_xyzs);

    // Free memory allocated for merging plan
    free((void *) obj->groups);
    free((void *) obj->polarities);

    // Free memory allocated for TDoAs
    free((void *) obj->tdoas);
    free((void *) obj->ranges);
    free((void *) obj->lookups);

    // Free memory allocated for cross-correlation
    free((void *) obj->cross_correlation);

    // Free memory allocated for FFT setup
    for (group_index = 0; group_index < obj->groups_count; group_index++) {
        free((void *) obj->frames_complex[group_index]);
        free((void *) obj->frames_real[group_index]);
        fftwf_destroy_plan(obj->iffts[group_index]);
    }
    free((void *) obj->frames_complex);
    free((void *) obj->frames_real);
    free((void *) obj->iffts);
    fftwf_cleanup();    

    // Free memory allocated for object
    free((void *) obj);

}

int smp_call(smp_obj * obj, const covs_obj * covs, ldoas_obj * ldoas) {

    unsigned int group_index;
    unsigned int pair_index;
    unsigned int bin_index;
    unsigned int pdoa_index;
    unsigned int offset;
    float energy;
    float energy_max;
    unsigned int pdoa_index_star;
    unsigned int lookup_index;    

    // Generate cross-correlations
    for (group_index = 0; group_index < obj->groups_count; group_index++) {

        // Reset all values to zero (this also handle interpolation)
        memset(obj->frames_complex[group_index], 0x00, sizeof(fftwf_complex) * ((obj->frame_size * obj->interpolation_rate)/2+1));

    }

    // For each pair in group, add cross-spectrum according to polarity
    for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {

        group_index = obj->groups[pair_index];

        if (obj->polarities[pair_index] == +1) {
            for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
                obj->frames_complex[group_index][bin_index][0] += covs->samples[pair_index][bin_index*2+0];
                obj->frames_complex[group_index][bin_index][1] += covs->samples[pair_index][bin_index*2+1];
            }
        }

        if (obj->polarities[pair_index] == -1) {
            for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
                obj->frames_complex[group_index][bin_index][0] += covs->samples[pair_index][bin_index*2+0];
                obj->frames_complex[group_index][bin_index][1] -= covs->samples[pair_index][bin_index*2+1];
            }
        }

    }

    // Process each group
    offset = 0;

    for (group_index = 0; group_index < obj->groups_count; group_index++) {

        // Perform IFFT for this group
        fftwf_execute(obj->iffts[group_index]);

        // Copy to supervector
        memcpy( &(obj->cross_correlation[offset+(obj->ranges[group_index]-1)/2]), 
                &(obj->frames_real[group_index][0]), 
                ((obj->ranges[group_index]+1)/2) * sizeof(float));

        memcpy( &(obj->cross_correlation[offset]), 
                &(obj->frames_real[group_index][obj->frame_size*obj->interpolation_rate-(obj->ranges[group_index]-1)/2]), 
                ((obj->ranges[group_index]-1)/2) * sizeof(float));

        // Increase offset pointer
        offset += obj->ranges[group_index];

    }

    // Scan the points
    energy_max = 0.0;
    lookup_index = 0;
    pdoa_index_star = 0;
    for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {

        energy = 0.0;
        for (group_index = 0; group_index < obj->groups_count; group_index++) {
            energy += obj->cross_correlation[obj->lookups[lookup_index]];
            lookup_index++;
        }

        if (energy > energy_max) {
            energy_max = energy;
            pdoa_index_star = pdoa_index;
        }

    }

    // Find corresponding point
    ldoas->xyzs[0] = obj->pdoas_xyzs[pdoa_index_star * 3 + 0];
    ldoas->xyzs[1] = obj->pdoas_xyzs[pdoa_index_star * 3 + 1];
    ldoas->xyzs[2] = obj->pdoas_xyzs[pdoa_index_star * 3 + 2];
    ldoas->es[0] = energy_max;

    return 0;

}

srp_obj * srp_construct(const unsigned int frame_size, const unsigned int interpolation_rate, const float sound_speed, const unsigned int sample_rate, const mics_obj * mics, const pdoas_obj * pdoas) {

    srp_obj * obj;

    unsigned int mic_index1;
    unsigned int mic_index2;
    unsigned int pair_index;
    unsigned int pdoa_index;
    float tdoa;
    unsigned int offset;

    // Allocate memory for the object
    obj = (srp_obj *) malloc(sizeof(smp_obj));

    // Copy relevant parameters
    obj->frame_size = frame_size;
    obj->interpolation_rate = interpolation_rate;
    obj->sound_speed = sound_speed;
    obj->sample_rate = sample_rate;
    obj->pdoas_count = pdoas->points_count;
    obj->pdoas_xyzs = (float *) malloc(sizeof(float) * pdoas->points_count * 3);
    memcpy(obj->pdoas_xyzs, pdoas->xyzs, sizeof(float) * pdoas->points_count * 3);
    
    // Compute vectors for each pair
    obj->pairs_count = (mics->channels_count * (mics->channels_count-1) / 2);
    obj->pairs_xyzs = (float *) malloc(sizeof(float) * obj->pairs_count * 3);
    pair_index = 0;
    for (mic_index1 = 0; mic_index1 < mics->channels_count; mic_index1++) {
        for (mic_index2 = (mic_index1+1); mic_index2 < mics->channels_count; mic_index2++) {
            obj->pairs_xyzs[pair_index*3+0] = mics->xyzs[mic_index1*3+0] - mics->xyzs[mic_index2*3+0];
            obj->pairs_xyzs[pair_index*3+1] = mics->xyzs[mic_index1*3+1] - mics->xyzs[mic_index2*3+1];
            obj->pairs_xyzs[pair_index*3+2] = mics->xyzs[mic_index1*3+2] - mics->xyzs[mic_index2*3+2];
            pair_index++;
        }
    }


    // Compute the TDoAs for the reference pair in each group
    obj->tdoas = (float *) malloc(sizeof(float) * obj->pairs_count * obj->pdoas_count);

    for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {        

        for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {

            // tdoa = -(fs/c) * (d.u)
            tdoa = 0.0;
            tdoa += obj->pdoas_xyzs[pdoa_index*3+0] * obj->pairs_xyzs[pair_index*3+0];
            tdoa += obj->pdoas_xyzs[pdoa_index*3+1] * obj->pairs_xyzs[pair_index*3+1];
            tdoa += obj->pdoas_xyzs[pdoa_index*3+2] * obj->pairs_xyzs[pair_index*3+2];
            tdoa *= (float) sample_rate;
            tdoa /= sound_speed;
            tdoa *= -1.0;

            obj->tdoas[pdoa_index * obj->pairs_count + pair_index] = tdoa;

        }

    }

    // Compute the range for each pair
    obj->ranges = (unsigned int *) malloc(sizeof(unsigned int) * obj->pairs_count); 
    for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {
        tdoa = 0.0;
        for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {
            if (fabs(obj->tdoas[pdoa_index * obj->pairs_count + pair_index]) > tdoa) {
                tdoa = fabs(obj->tdoas[pdoa_index * obj->pairs_count + pair_index]);
            }
        }
        obj->ranges[pair_index] = ((unsigned int) ceilf(tdoa * (float) interpolation_rate)) * 2 + 1;
    }

    // Compute number of samples in the super cross-correlation vector
    offset = 0;
    for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {
        offset += obj->ranges[pair_index];
    }

    // Generate the lookup table
    offset = 0;
    obj->lookups = (unsigned int *) malloc(sizeof(unsigned int) * obj->pairs_count * obj->pdoas_count);    
    for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {
        for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) { 
            tdoa = obj->tdoas[pdoa_index * obj->pairs_count + pair_index];
            obj->lookups[pdoa_index * obj->pairs_count + pair_index] = offset + (obj->ranges[pair_index]-1)/2 + ((unsigned int) roundf(tdoa * (float) interpolation_rate));
        }
        offset += obj->ranges[pair_index];
    }

    // Create a cross-correlation supervector
    obj->cross_correlation = (float *) malloc(sizeof(float) * offset);

    // Create the FFT setup
    obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * ((frame_size * interpolation_rate)/2+1));
    obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size * interpolation_rate);
    obj->ifft = fftwf_plan_dft_c2r_1d(frame_size * interpolation_rate, obj->frame_complex, obj->frame_real, FFTW_ESTIMATE);    

    return obj;

}

void srp_destroy(srp_obj * obj) {

    // Free memory allocated for doas
    free((void *) obj->pdoas_xyzs);

    // Free memory allocated to store pairs differences
    free((void *) obj->pairs_xyzs);

    // Free memory allocated for TDoAs
    free((void *) obj->tdoas);
    free((void *) obj->ranges);
    free((void *) obj->lookups);

    // Free memory allocated for cross-correlation
    free((void *) obj->cross_correlation);

    // Free memory allocated for FFT setup
    free((void *) obj->frame_complex);
    free((void *) obj->frame_real);
    fftwf_destroy_plan(obj->ifft);
    fftwf_cleanup();    

    // Free memory allocated for object
    free((void *) obj);

}

int srp_call(srp_obj * obj, const covs_obj * covs, ldoas_obj * ldoas) {

    unsigned int pair_index;
    unsigned int bin_index;
    unsigned int pdoa_index;
    unsigned int offset;
    float energy;
    float energy_max;
    unsigned int pdoa_index_star;
    unsigned int lookup_index;

    offset = 0;
    
    // Generate cross-correlations
    for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {

        // Reset all values to zero (this also handle interpolation)
        memset(obj->frame_complex, 0x00, sizeof(fftwf_complex) * ((obj->frame_size * obj->interpolation_rate)/2+1));

        // Copy the input frequency signal
        for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
            obj->frame_complex[bin_index][0] = covs->samples[pair_index][bin_index*2+0];
            obj->frame_complex[bin_index][1] = covs->samples[pair_index][bin_index*2+1];
        }
                
        // Perform IFFT for this group
        fftwf_execute(obj->ifft);

        // Copy to supervector
        memcpy( &(obj->cross_correlation[offset+(obj->ranges[pair_index]-1)/2]), 
                &(obj->frame_real[0]), 
                ((obj->ranges[pair_index]+1)/2) * sizeof(float));

        memcpy( &(obj->cross_correlation[offset]), 
                &(obj->frame_real[obj->frame_size*obj->interpolation_rate-(obj->ranges[pair_index]-1)/2]), 
                ((obj->ranges[pair_index]-1)/2) * sizeof(float));

        // Increase offset pointer
        offset += obj->ranges[pair_index];

    }


    // Scan the points
    energy_max = 0.0;
    lookup_index = 0;
    pdoa_index_star = 0;
    for (pdoa_index = 0; pdoa_index < obj->pdoas_count; pdoa_index++) {

        energy = 0.0;
        for (pair_index = 0; pair_index < obj->pairs_count; pair_index++) {
            energy += obj->cross_correlation[obj->lookups[lookup_index]];
            lookup_index++;
        }

        if (energy > energy_max) {
            energy_max = energy;
            pdoa_index_star = pdoa_index;
        }

    }

    // Find corresponding point
    ldoas->xyzs[0] = obj->pdoas_xyzs[pdoa_index_star * 3 + 0];
    ldoas->xyzs[1] = obj->pdoas_xyzs[pdoa_index_star * 3 + 1];
    ldoas->xyzs[2] = obj->pdoas_xyzs[pdoa_index_star * 3 + 2];
    ldoas->es[0] = energy_max;

    return 0;

}

csv_obj * csv_construct(const char * file_name, const unsigned int channels_count) {

    csv_obj * obj;

    // Allocate memory for object
    obj = (csv_obj *) malloc(sizeof(csv_obj));

    // Create CSV file
    obj->file_pointer = fopen(file_name, "w");

    // Copy relevant parameters    
    obj->frame_index = 0;
    obj->channels_count = channels_count;

    // Write CSV header to the file
    fprintf(obj->file_pointer, "frame_index, x, y, z\n");

    // Return pointer to object
    return obj;

}

void csv_destroy(csv_obj * obj) {

    // Close CSV file
    fclose(obj->file_pointer);

    // Deallocate memory for object
    free((void *) obj);

}

int csv_write(csv_obj * obj, ldoas_obj * ldoas) {

    obj->frame_index++;

    // Write a new line with the frame index first
    fprintf(obj->file_pointer, "%u, %+1.4f, %+1.4f, %+1.4f\n", obj->frame_index, ldoas->xyzs[0], ldoas->xyzs[1], ldoas->xyzs[2]);

    // By default return 0 when returns without error        
    return 0;

}