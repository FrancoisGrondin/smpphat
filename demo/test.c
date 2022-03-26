#include <smpphat/signal.h>
#include <smpphat/system.h>

#include <time.h>

int main(int argc, char * argv[]) {

	mics_obj * mics;
	pdoas_obj * pdoas;
	covs_obj * covs;
	ldoas_obj * ldoas;
	smp_obj * smp;
	srp_obj * srp;

	clock_t begin;
	clock_t end;
	float time_spent;
	float srp_time;
	float smp_time;
	float total_time;

	const unsigned int levels_count = 4;
	const unsigned int frame_size = 512;
	const unsigned int interpolation_rate = 4;
	const float sound_speed = 343.0;
	const unsigned int sample_rate = 16000;
	const unsigned int ldoas_count = 1;

	const unsigned int mic_arrays_count = 4;
	const char mic_arrays[][32] = { "respeaker_usb", "respeaker_core", "minidsp_uma", "matrix_creator" };

	unsigned int mic_arrays_index;
	unsigned int iteration_index;
	unsigned int iterations_count;

	for (mic_arrays_index = 0; mic_arrays_index < mic_arrays_count; mic_arrays_index++) {

		iterations_count = 1;

		mics = mics_construct(mic_arrays[mic_arrays_index]);
	    pdoas = pdoas_construct(levels_count, "half");
	    covs = covs_construct(mics->channels_count, frame_size);
	    ldoas = ldoas_construct(ldoas_count);
		smp = smp_construct(frame_size, interpolation_rate, sound_speed, sample_rate, mics, pdoas);
		srp = srp_construct(frame_size, interpolation_rate, sound_speed, sample_rate, mics, pdoas);

		while(1) {

			total_time = 0.0;

			begin = clock();
			for (iteration_index = 0; iteration_index < iterations_count; iteration_index++) {
				srp_call(srp, covs, ldoas);		
			}
			end = clock();
			time_spent = (float) (end - begin) / CLOCKS_PER_SEC;
			srp_time = 1000.0 * time_spent / ((float) iterations_count);

			total_time += time_spent;

			begin = clock();
			for (iteration_index = 0; iteration_index < iterations_count; iteration_index++) {
				smp_call(smp, covs, ldoas);		
			}
			end = clock();
			time_spent = (float) (end - begin) / CLOCKS_PER_SEC;
			smp_time = 1000.0 * time_spent / ((float) iterations_count);

			total_time += time_spent;

			if (total_time > 2.0) {
				break;
			}
			else {
				iterations_count *= 2;
			}

		}

		mics_destroy(mics);
		pdoas_destroy(pdoas);
		covs_destroy(covs);
		ldoas_destroy(ldoas);
		smp_destroy(smp);
		srp_destroy(srp);

		printf("%s\n", mic_arrays[mic_arrays_index]);
		printf("Srp: %f msec/frame\n", srp_time);
		printf("Smp: %f msec/frame\n", smp_time);
		printf("\n");

	}

	return 0;

}