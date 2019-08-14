// fast and slow convolution routines, with autoscaling  + lamp shift -- header file
// Version under construction!!!


#ifndef N_CHANNELS
 #define N_CHANNELS 8192
#endif

#ifndef DELTA_T
 #define DELTA_T 0.05
#endif

/* rescaling old */
void rescale(double*, double*, double*, int, int);

/* rescaling new */
void rescale_w(double*, double*, double*, double*, int, int);

/* rescaling new + bg */
void rescale_w_bg(double*, double*, double*, double, double*, int, int);

/* fast convolution */
void fconv(double*, double*, double*, int, int, int);

/* fast convolution, high repetition rate */
void fconv_per(double*, double*, double*, int, int, int, int, double);

/* fast convolution, high repetition rate, with convolution stop */
void fconv_per_cs(double*, double*, double*, int, int, int, double, int);
		   
/* fast convolution with reference compound */
void fconv_ref(double*, double*, double*, int, int, int, double);

/* slow convolution */
void sconv(double*, double*, double*, int, int);

/* shifting lamp */
void shift_lamp(double*, double*, double, int);
