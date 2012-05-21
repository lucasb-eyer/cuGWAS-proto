#ifndef FGLS_CHOL_GPU_H
#define FGLS_CHOL_GPU_H

int fgls_chol_gpu(
        int n, int p, int m, int t, int wXL, int wXR,
        int x_b, int y_b, int num_threads, // y_b not used
        char *Phi_path, char *h2_path, char *sigma2_path,
        char *XL_path, char *XR_path, char *Y_path,
        char *B_path, char *V_path
);

#endif // FGLS_CHOL_GPU_H
