/* Calibration table for AYB, produced by ayb-recal
 * Bpert cycle, end2 (sampled)
 * Created: Tue Oct 25 15:48:52 2011
 */

#ifndef HAS_CALIBRATION_TABLE
#define HAS_CALIBRATION_TABLE

real_t calibration_intercept =  14.395204;

real_t calibration_scale     =  0.319750;

real_t calibration_baseprior_adj[] = {
 1.666828,  0.124661, -0.736526,  0.753566,	// A?
-0.756845,  0.229444,  1.342972, -1.412537,	// C?
-1.913303, -2.809816, -1.074479, -0.092334,	// G?
-1.306052,  1.913722,  2.128846,  1.941854,	// T?
};

real_t calibration_basenext_adj[] = {
-0.415869,  0.994301,  0.542031, -0.731049,	// ?A
-4.266027,  1.571106,  4.803866, -3.606623,	// ?C
-4.171984,  1.858490,  3.045182, -2.334754,	// ?G
-2.906095,  0.984550,  2.831418,  1.801458,	// ?T
};

real_t calibration_priorbasenext_adj[] = {
 2.309614, -1.819272, -1.924556,  1.289015,	// A?A
 0.409760, -0.109266, -0.570290,  0.282671,	// C?A
 2.726730, -0.332101, -0.145963,  2.156758,	// G?A
 0.437819, -0.450728,  0.119350,  0.861007,	// T?A
 0.634546, -0.432383, -0.749551,  0.473195,	// A?C
 0.681140, -5.466834, -2.867783,  4.311379,	// C?C
 2.401949,  1.101690, -2.655685,  4.021926,	// G?C
 3.506246, -2.957882, -8.706013,  4.380209,	// T?C
 2.259817, -0.420950,  0.125216,  1.637828,	// A?G
 3.717892, -3.810383, -4.636865,  3.910522,	// C?G
 1.436970,  2.027164, -4.305729,  0.329824,	// G?G
 1.460135, -3.628185, -2.905160,  1.269832,	// T?G
 0.679946, -0.038763,  0.027432,  1.189412,	// A?T
 3.476472, -2.759751, -2.765033,  3.672195,	// C?T
-0.060261,  2.632756, -1.484128, -0.591495,	// G?T
 2.233145, -2.970263, -5.150166,  0.523869,	// T?T
};

#endif  /*  HAS_CALIBRATION_TABLE  */
