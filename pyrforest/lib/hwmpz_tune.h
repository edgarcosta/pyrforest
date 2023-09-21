// File generated by profile_hwmpz
#ifndef _HWMPZ_TUNE_INCLUDE_
#define _HWMPZ_TUNE_INCLUDE_

#ifdef __cplusplus
extern "C" {
#endif
#define HWMPZ_TUNE_MAX_D	40	// for performance reasons this bound is NOT verified, caller must verify

extern int mpz_row_fft_crossovers[HWMPZ_TUNE_MAX_D+1];
static inline long mpz_row_fft_crossover (int d) { return mpz_row_fft_crossovers[d <=  HWMPZ_TUNE_MAX_D ? d : HWMPZ_TUNE_MAX_D]; }

extern int mpz_mat_fft_crossovers[HWMPZ_TUNE_MAX_D+1];
static inline long mpz_mat_fft_crossover (int d) { return mpz_mat_fft_crossovers[d <=  HWMPZ_TUNE_MAX_D ? d : HWMPZ_TUNE_MAX_D]; }

extern int mpz_mod_fft_crossovers[HWMPZ_TUNE_MAX_D*HWMPZ_TUNE_MAX_D+1];
static inline long mpz_mod_fft_crossover (int d) { return mpz_mod_fft_crossovers[d <=  HWMPZ_TUNE_MAX_D ? d : HWMPZ_TUNE_MAX_D]; }

#endif
