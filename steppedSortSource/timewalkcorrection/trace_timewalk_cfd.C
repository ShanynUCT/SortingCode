#include <stdio.h>
#include <stdlib.h>


// To compile: gcc -o cfd_program trace_timewalk_cfd.c 
// To run: ./cfd_program <input_file>.dat

#define TRACE_SIZE 3999 // Adjust as needed for your file
#define FIRST_TRACES 400 // Number of traces to use for subtraction

long double CFD(int T[TRACE_SIZE], double S, int D);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    int trace[TRACE_SIZE];
    double scale = 0.01; // Example scale
    int delay = 1;       // Example delay

    // Open the .dat file
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return 1;
    }

    // Read the file line by line
    int index, value;
    for (int i = 0; i < TRACE_SIZE; i++) {
        if (fscanf(file, "%d %d", &index, &value) != 2) {
            fprintf(stderr, "Error reading file %s at line %d\n", filename, i + 1);
            fclose(file);
            return 1;
        }
        trace[i] = value;
    }

    // Close the file
    fclose(file);

    // Calculate the average amplitude of the first 400 traces
    int sum = 0;
    for (int i = 0; i < FIRST_TRACES; i++) {
        sum += trace[i];
    }
    double average = (double)sum / FIRST_TRACES;

    // Subtract the average amplitude from the whole dataset
    for (int i = 0; i < TRACE_SIZE; i++) {
        trace[i] -= average;
        trace[i] *= -1;
    }

    // Output data for traces 400-800 to a file
    FILE *output_file = fopen("output.dat", "w");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening output file\n");
        return 1;
    }
    for (int i = 400; i < 800; i++) {
        fprintf(output_file, "%d %d\n", i, trace[i]);
    }
    fclose(output_file);

    return 1;
}

long double CFD(int T[TRACE_SIZE], double S, int D) {
    double T_sample = 2;
    double T_coarse = 0;
    double C_a, C_b;
    long double ZP;
    for (int j = D; j < TRACE_SIZE; j++) {
        C_a = (double)(S * T[j] - T[j - D] + (1 - S) * T[0]);
        if (j == D) C_b = C_a;
        if (C_b > 0. && C_a < 0.) {
            ZP = T_coarse + T_sample * (C_b / (C_b - C_a));
            C_b = C_a;
            printf("ZP: %Lf\n", ZP);
            return ZP;
        } else if (C_b > 0. && C_a > 0.) {
            T_coarse += T_sample;
            C_b = C_a;
        }
    }
    return -1;
}
