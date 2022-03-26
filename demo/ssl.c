#include <smpphat/signal.h>
#include <smpphat/system.h>

#include <unistd.h>
#include <string.h>
#include <math.h>

int main(int argc, char * argv[]) {

    int opt;
    char method;
    char wav_filename[1024];
    char csv_filename[1024];
    char mic_filename[1024];

    wav_obj * wav;
    stft_obj * stft;
    scmphat_obj * scmphat;
	smp_obj * smp;
	srp_obj * srp;
	csv_obj * csv;

    hops_obj * hops;
    freqs_obj * freqs;
    covs_obj * covs;
	ldoas_obj * ldoas;    

	mics_obj * mics;
	pdoas_obj * pdoas;

    const unsigned int hop_size = 256;
    const unsigned int frame_size = 512;
    const unsigned int sample_rate = 16000;
    const unsigned int interpolation_rate = 4;
    const float sound_speed = 343.0;
    const float alpha = 1.0;
    const unsigned int levels_count = 4;

    method = 0x00;

    while((opt = getopt(argc, argv, "a:i:m:o:")) != -1) 
    { 
        switch(opt) 
        { 
            
        	case 'a':

        		if (strcmp(optarg, "respeaker_usb") == 0) {
        			strcpy(mic_filename, "respeaker_usb");
        		}
        		if (strcmp(optarg, "respeaker_core") == 0) {
        			strcpy(mic_filename, "respeaker_core");
        		}
        		if (strcmp(optarg, "minidsp_uma") == 0) {
        			strcpy(mic_filename, "minidsp_uma");
        		}
        		if (strcmp(optarg, "matrix_creator") == 0) {
        			strcpy(mic_filename, "matrix_creator");
        		}

        		break;

            case 'i':

                strcpy(wav_filename, optarg);
                break;

            case 'm':

            	if (strcmp(optarg, "srp") == 0) {
            		method = 'r';
            	}
            	if (strcmp(optarg, "smp") == 0) {
            		method = 'm';
            	}
            	break;

            case 'o':

                strcpy(csv_filename, optarg);
                break;

        } 
    } 

    if( access(wav_filename, F_OK) != 0 ) {
        printf("File %s does not exists\n", wav_filename);
        exit(-1);
    }

    wav = wav_construct(wav_filename, hop_size);

    if (wav->sample_rate != sample_rate) {
        printf("Invalid sample rate (was expecting %u, got %u).\n", sample_rate, wav->sample_rate);
        exit(-1);
    }

    mics = mics_construct(mic_filename);
    pdoas = pdoas_construct(levels_count, "half");

    if (wav->num_channels != mics->channels_count) {
    	printf("Mismatch between number of channels in wav file (%u) and number of microphones in array %s (%u).\n", wav->num_channels, mic_filename, mics->channels_count);
    	exit(-1);
    }

    stft = stft_construct(wav->num_channels, frame_size, hop_size);
    scmphat = scmphat_construct(wav->num_channels, frame_size, alpha);
	smp = smp_construct(frame_size, interpolation_rate, sound_speed, sample_rate, mics, pdoas);
	srp = srp_construct(frame_size, interpolation_rate, sound_speed, sample_rate, mics, pdoas);
	csv = csv_construct(csv_filename, wav->num_channels);

    hops = hops_construct(wav->num_channels, hop_size);
    freqs = freqs_construct(wav->num_channels, frame_size);
    covs = covs_construct(wav->num_channels, frame_size);
    ldoas = ldoas_construct(1);

    while (wav_read(wav, hops) == 0) {

        stft_call(stft, hops, freqs);
        scmphat_call(scmphat, freqs, covs);

        if (method == 'r') {
			srp_call(srp, covs, ldoas);    	
        }
        if (method == 'm') {
        	smp_call(smp, covs, ldoas);
        }

	    csv_write(csv, ldoas);

    }

    stft_destroy(stft);
    scmphat_destroy(scmphat);
    smp_destroy(smp);
    srp_destroy(srp);
    csv_destroy(csv);

    hops_destroy(hops);
    freqs_destroy(freqs);
    covs_destroy(covs);
    ldoas_destroy(ldoas);

	return 0;

}