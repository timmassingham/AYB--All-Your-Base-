/* Calibration table for AYB, produced by ayb-recal
 * NA19240_BGI end2
 * Created: Thu Sep  1 09:43:15 2011
 */

#ifndef HAS_CALIBRATION_TABLE
#define HAS_CALIBRATION_TABLE

real_t calibration_intercept =  12.338324;

real_t calibration_scale     =  0.327644;

real_t calibration_baseprior_adj[] = {
 1.423309, -0.598516,  0.004191,  0.284387,	// A?
 1.204328, -0.451003, -3.355595,  1.224770,	// C?
 1.325993, -3.240346, -0.107461,  1.296676,	// G?
 0.189029, -0.021458, -0.812684,  1.634378,	// T?
};

real_t calibration_basenext_adj[] = {
 1.503320, -0.752066, -0.125230,  0.009147,	// ?A
 0.851348, -0.069675, -3.487903,  1.842629,	// ?C
 1.936539, -3.352013,  0.317069,  0.904239,	// ?G
-0.148548, -0.137569, -0.975485,  1.684197,	// ?T
};

real_t calibration_priorbasenext_adj[] = {
 0.298461,  0.961372,  0.657088, -0.383188,	// A?A
-1.594070,  1.773384,  1.768246, -1.689642,	// C?A
-1.779165,  0.825758,  0.942141, -1.487750,	// G?A
-1.067884,  0.750807,  0.904073, -0.879631,	// T?A
-1.761967,  1.660403,  1.766007, -1.860658,	// A?C
-0.800411, -0.356713,  0.041550, -0.228235,	// C?C
-0.349182,  2.085897,  1.695911, -0.455861,	// G?C
-1.231098,  0.921736,  0.768081, -1.895458,	// T?C
-1.482486,  0.941595,  0.925768, -1.445868,	// A?G
-0.499814,  2.021779,  1.789027, -0.356596,	// C?G
-0.408331, -0.479714,  0.288190, -0.769400,	// G?G
-1.752027,  1.827662,  1.268563, -1.868348,	// T?G
-1.196666,  0.747952,  0.922685, -0.750498,	// A?T
-1.248363,  0.872872,  0.672725, -2.165739,	// C?T
-1.605980,  1.879381,  1.345307, -1.727201,	// G?T
-0.091649,  0.811117,  1.330831,  0.203225,	// T?T
};

#endif  /*  HAS_CALIBRATION_TABLE  */