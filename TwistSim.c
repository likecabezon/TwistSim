#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define SPEED_OF_LIGHT 2.99792e8  // Speed of light in m/s
#define MU_0 1.256637e-6
#define EPS_0 8.8541878188e-12
// Uncomment if your math library doesnt include M_PI
//#define M_PI 3.14159265359

typedef struct {
    double x;
    double y;
    double z;
} Point3D;

Point3D* generate_coil_points(double radius, double length, int n_turns, double wire_rad, double min_freq, double max_freq, double extra_height, int *num_points,
                              double *segment_L, int segmentation_mode, int segmentation_data) {
  // Segmentation: 0 = auto, 1 = segments per turn, 2 = segment number,
  // 3 = segment wavelength density(number of segments per wavelength)
  // Constants
  const int auto_min_wavelength_density = 20, auto_segments_per_turn = 16;
  // Total angular displacement over the entire coil (in radians)
  double max_angle = n_turns * 2 * M_PI;

  // Length covered per unit of angular displacement
  double normal_length_per_angle = length / max_angle;

  double segment_length;

  if ((2.0 * wire_rad * n_turns) > length){
    printf("Warning: Wires overlapped by the excessive wire radius number of turns or not enough coil length\n");
  }


  if(segmentation_mode == 3){
    segment_length = SPEED_OF_LIGHT/(min_freq * segmentation_data);

    printf("mode 3 choosen\n");

    // Variables for Newton-Raphson iteration
    int max_iterations = 20;
    double previous_guess = 0.0;
    double actual_guess = M_PI / 2.0;  // Initial guess for angular step size
    int i = 0;
    double relative_epsilon = 1e-9;

    // Newton-Raphson iteration to find the correct step size (angle per segment)
    while (fabs((actual_guess - previous_guess) / actual_guess) > relative_epsilon && i < max_iterations) {
      previous_guess = actual_guess;

      // Newton-Raphson step
      actual_guess = actual_guess-(pow(sin(previous_guess / 2)*radius * 2, 2) + pow(normal_length_per_angle * actual_guess, 2) - segment_length * segment_length) / (2*(radius * radius * sin(actual_guess) + normal_length_per_angle * normal_length_per_angle * actual_guess));

      i++;
    }

  // If the iteration doesn't converge, print a warning
    if (i == max_iterations) {
      printf("Warning: Newton-Raphson did not converge within the maximum iterations\n");
    }

    // Calculate the number of steps (discretized segments)
    *num_points = (int)(max_angle / actual_guess) + 1;
  }

  if(segmentation_mode == 2){
    printf("mode 2 choosen\n");
    *num_points = segmentation_data +1;
  }

  if(segmentation_mode == 1){
    printf("mode 1 choosen\n");
    *num_points = n_turns * segmentation_data + 1;
  }

  if(segmentation_mode == 0){
    printf("mode 0 choosen\n");
    *num_points = n_turns * auto_segments_per_turn +1;
    segment_length = sqrt(pow(2*radius*sin(max_angle/(2*(*num_points -1))),2) + pow(normal_length_per_angle*max_angle/(*num_points -1),2));

    if(segment_length > SPEED_OF_LIGHT/(max_freq * auto_min_wavelength_density)){

      segment_length = SPEED_OF_LIGHT / (max_freq * auto_min_wavelength_density);

      // Variables for Newton-Raphson iteration
      int max_iterations = 20;
      double previous_guess = 0.0;
      double actual_guess = M_PI / 2.0;  // Initial guess for angular step size
      int i = 0;
      double relative_epsilon = 1e-9;

      // Newton-Raphson iteration to find the correct step size (angle per segment)
      while (fabs((actual_guess - previous_guess) / actual_guess) > relative_epsilon && i < max_iterations) {
        previous_guess = actual_guess;

        // Newton-Raphson step
        actual_guess = actual_guess-(pow(sin(previous_guess / 2)*radius * 2, 2) + pow(normal_length_per_angle * actual_guess, 2) - segment_length * segment_length) / (2*(radius * radius * sin(actual_guess) + normal_length_per_angle * normal_length_per_angle * actual_guess));

        i++;
      }

      // If the iteration doesn't converge, print a warning
      if (i == max_iterations) {
        printf("Warning: Newton-Raphson did not converge within the maximum iterations\n");
      }

      // Calculate the number of steps (discretized segments)
      *num_points = (int)(max_angle / actual_guess) + 2;
    }
  }

  // Aprox length of wire needed based on current segmentation
  *segment_L = sqrt(pow(2*radius*sin(max_angle/(2*(*num_points -1))),2) + pow(normal_length_per_angle*max_angle/(*num_points -1),2));

  double wire_length =  (*segment_L) * (*num_points -1) + 2*extra_height;

  printf("num of points = %d", *num_points);


  // Allocate memory for the points
  Point3D *points = (Point3D *)malloc((*num_points) * sizeof(Point3D));
  if (!points) {
    printf("Error: Memory allocation failed\n");
    return NULL;
  }

  // Generate theta values for angular positions along the coil
  double *theta_vec = (double *)malloc((*num_points) * sizeof(double));
  if (!theta_vec) {
    printf("Error: Memory allocation failed\n");
    free(points);
    return NULL;
  }

  if ((n_turns % 2) == 0) {
    for (int j = 0; j < *num_points; j++) {
      theta_vec[j] = -(max_angle + M_PI) / 2.0 + j * (max_angle) / (*num_points - 1);
    }
  }
  else {
      for (int j = 0; j < *num_points; j++) {
          theta_vec[j] = (M_PI - max_angle) / 2.0 + j * (max_angle) / (*num_points - 1);
      }
  }

  // Generate the corresponding z positions (normal length along the coil)
  double *normal_len_vec = (double *)malloc((*num_points) * sizeof(double));
  if (!normal_len_vec) {
      printf("Error: Memory allocation failed\n");
      free(points);
      free(theta_vec);
      return NULL;
  }

  for (int j = 0; j < *num_points; j++) {
      normal_len_vec[j] = -length / 2.0 + j * length / (*num_points-1);
  }

  // Fill the point array with x, y, z coordinates
  for (int j = 0; j < *num_points; j++) {
      points[j].x = normal_len_vec[j];
      points[j].y = radius * cos(theta_vec[j]);
      points[j].z = radius * (sin(theta_vec[j]) + 1) + extra_height;
  }

  // Free the temporary arrays-
  free(theta_vec);
  free(normal_len_vec);

  // Return the array of points
  return points;
}


// cambiar las funciones de resistencia de los segmentos con las ecuaciones de skin y proximity
int create_nec_file(const char* filename, Point3D* points, int num_points, double wire_radius, double freq, double conductivity,
                    double complex **currents, double (*mag_interactions[])[8], double proximity_effect_constant,
                    double proximity_effect_constant_bases, int compute_effects, int *num_iterations, int max_iterations,
                    double conver_max_rel_err) {
  /*    Compute_effects modes
      0 No resistance in wire                                   -Default
      1 Wire conductivity                                       -No computation
      2 Skin effect equivalent impedance                        -Direct solution
      3 Skin effect + proximity effect equivalent impedance     -Needs iterations(this function becomes part of a loop)
  */


  FILE* fp = fopen(filename, "w");
  if (fp == NULL) {
    perror("Error opening file");
    return -1;
  }

  // NEC2 comment lines
  fprintf(fp, "CM Coil over Ground Plane\n");
  fprintf(fp, "CE\n");



  // Add the excitation element from the ground plane (z = 0) to the first point
  fprintf(fp, "GW 1 1 %lf %lf %lf %lf %lf %lf %lf\n",
          points[0].x, points[0].y, 0.0,  // Last point of coil
          points[0].x, points[0].y, points[0].z,  // Connect to ground plane (z = 0)
          wire_radius);

  // Define the coil geometry using GW cards
  for (int i = 2; i <= num_points; i++) {
    fprintf(fp, "GW %d 1 %lf %lf %lf %lf %lf %lf %lf\n",
            i,  // Wire tag number (incremental for each wire)
            points[i-2].x, points[i-2].y, points[i-2].z,  // Starting point (x1, y1, z1)
            points[i - 1].x, points[i - 1].y, points[i - 1].z,  // Ending point (x2, y2, z2)
            wire_radius);  // Wire radius
  }

  // Add the short from the last point to the ground plane (z = 0)
  fprintf(fp, "GW %d 1 %lf %lf %lf %lf %lf %lf %lf\n",
          num_points + 1,
          points[num_points - 1].x, points[num_points - 1].y, points[num_points - 1].z,  // Last point of coil
          points[num_points - 1].x, points[num_points - 1].y, 0.0,  // Connect to ground plane (z = 0)
          wire_radius);

  // End of geometry definition
  fprintf(fp, "GE 1\n");


  if(compute_effects == 1){
    // Add conductivity to the wires
    for (int i = 1; i < num_points + 2; i++){
      fprintf(fp, "LD 5 0 %d 0 %lf\n", i, conductivity);
    }
  }
  else if(compute_effects == 2){
    double skin_depth = sqrt((sqrt(1 + pow(2*M_PI*freq*EPS_0/conductivity,2)) + 2*M_PI*freq*EPS_0/conductivity)/(freq*M_PI*MU_0*conductivity));
    double loss_coeff = (2*wire_radius - skin_depth*(1 - exp(-2*wire_radius/conductivity)))/(wire_radius*(1 - 2*exp(-wire_radius/skin_depth)*cos(wire_radius/skin_depth) + exp(-2*wire_radius/skin_depth)));
    double resistance = loss_coeff * points[0].z/conductivity;
    fprintf(fp, "LD 4 0 1 0 %lf 0\n", resistance);

    resistance = loss_coeff * sqrt(pow(points[0].x - points[1].x,2) + pow(points[0].y - points[1].y,2) + pow(points[0].z - points[1].z,2))/conductivity;
    for(int i = 2; i < num_points + 1; i++){
      fprintf(fp, "LD 4 0 %d 0 %lf 0\n", i, resistance);
    }

    resistance = loss_coeff * points[num_points - 1].z/conductivity;
    fprintf(fp, "LD 4 0 %d 0 %lf 0\n", num_points + 1, resistance);
  }
  else if(compute_effects == 3){
    double len_segs = sqrt(pow(points[1].x - points[0].x,2) + pow(points[1].y - points[0].y,2) + pow(points[1].z - points[0].z,2));
    double len_bases = points[0].z;

    if((*num_iterations) == 0){
        //completar para la primera iteracion(corrientes iguales con efecto skin)
        double skin_depth = sqrt((sqrt(1 + pow(2*M_PI*freq*EPS_0/conductivity,2)) + 2*M_PI*freq*EPS_0/conductivity)/(freq*M_PI*MU_0*conductivity));
        double skin_loss_coeff = (2*wire_radius - skin_depth*(1 - exp(-2*wire_radius/conductivity)))/(wire_radius*(1 - 2*exp(-wire_radius/skin_depth)*cos(wire_radius/skin_depth) + exp(-2*wire_radius/skin_depth)));
        double skin_curr = conductivity/(skin_loss_coeff*((num_points-1)*len_segs + 2*len_bases));
        double cum_resistance = 0;

        for(int i = 0; i <= num_points; i++){
            double len;
            double complex b_field_x = 0, b_field_y = 0, b_field_z = 0;
            double b_tot_squared, prox_current, eq_resistance;
            for(int j = 0; j <= num_points; j++){
                if(i != j){
                    b_field_x += skin_curr * (mag_interactions[i][j][0]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][3]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));


                    b_field_y += skin_curr * (mag_interactions[i][j][1]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][4]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));

                    b_field_z += skin_curr * (mag_interactions[i][j][2]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][5]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));
                }
                else{
                    b_field_x += skin_curr * mag_interactions[i][j][3]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);

                    b_field_y += skin_curr * mag_interactions[i][j][4]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);

                    b_field_z += skin_curr * mag_interactions[i][j][5]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);
                }

            }
            //printf("\n%e %e %e\n", cabs(b_field_x),cabs(b_field_y),cabs(b_field_z));
            b_tot_squared = pow(creal(b_field_x),2) + pow(creal(b_field_y),2) + pow(creal(b_field_z),2) + pow(cimag(b_field_x),2) + pow(cimag(b_field_y),2) + pow(cimag(b_field_z),2);

            prox_current = M_PI*b_tot_squared*pow(freq*M_PI*conductivity*wire_radius*wire_radius,2);

            if((i == 0) || (i == num_points)){
                len = len_bases;
            }
            else{
                len = len_segs;
            }

            eq_resistance = len * (skin_loss_coeff + prox_current/cpow(skin_curr,2))/conductivity;
            fprintf(fp, "LD 4 0 %d 0 %lf 0\n", i+1, eq_resistance);
            cum_resistance += eq_resistance;
            printf("\n%e",eq_resistance);
        }
        double complex curr_stimated = 1/cum_resistance;
        for(int i = 0; i <= num_points; i++){
            currents[0][i] = curr_stimated;
            currents[1][i] = 0;
        }
        (*num_iterations)++;
    }
    else{
      double rel_error = 0;
      double mean_squared_coil_curr = 0;
      //printf("\n%d",*num_iterations);
      for (int i = 0; i <= num_points; i++) {
        mean_squared_coil_curr += pow(cabs(currents[0][i]),2)/(num_points + 1);
        // Compute relative squared error for this segment
        if (cabs(currents[0][i]) > 1e-7) { // Avoid division by zero

          rel_error += cabs(currents[0][i] - currents[1][i])/cabs(currents[0][i]);
          //printf("Rel_error = %e; currents[0][%d] = %e + j%e; currents[1][%d] = %e + j%e\n", rel_error, i, creal(currents[0][i]), cimag(currents[0][i]), i, creal(currents[1][i]), cimag(currents[1][i]));
        }
        else if(cabs(currents[1][i]) > 1e-7){ // Avoid division by zero

          rel_error += cabs(currents[0][i] - currents[1][i])/cabs(currents[1][i]);

          //printf("Rel_error = %e; currents[0][%d] = %e + j%e; currents[1][%d] = %e + j%e\n", rel_error, i, creal(currents[0][i]), cimag(currents[0][i]), i, creal(currents[1][i]), cimag(currents[1][i]));
        }
    }

      if(((*num_iterations) < max_iterations) && (rel_error >= conver_max_rel_err)){

        double skin_depth = sqrt((sqrt(1 + pow(2*M_PI*freq*EPS_0/conductivity,2)) + 2*M_PI*freq*EPS_0/conductivity)/(freq*M_PI*MU_0*conductivity));
        double skin_loss_coeff = (2*wire_radius - skin_depth*(1 - exp(-2*wire_radius/conductivity)))/(wire_radius*(1 - 2*exp(-wire_radius/skin_depth)*cos(wire_radius/skin_depth) + exp(-2*wire_radius/skin_depth)));

        for(int i = 0; i <= num_points; i++){
            double len;
            double complex b_field_x = 0, b_field_y = 0, b_field_z = 0;
            double b_tot_squared, prox_current, eq_resistance;
            for(int j = 0; j <= num_points; j++){
                if(i != j){
                    b_field_x += (1*currents[0][j] + 0*currents[1][j]) * (mag_interactions[i][j][0]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][3]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));


                    b_field_y += (1*currents[0][j] + 0*currents[1][j]) * (mag_interactions[i][j][1]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][4]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));

                    b_field_z += (1*currents[0][j] + 0*currents[1][j]) * (mag_interactions[i][j][2]*cexp(mag_interactions[i][j][6]*2*M_PI*freq*I/SPEED_OF_LIGHT)
                            + mag_interactions[i][j][5]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT));
                }
                else{
                    b_field_x += (1*currents[0][j] + 0*currents[1][j]) * mag_interactions[i][j][3]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);

                    b_field_y += (1*currents[0][j] + 0*currents[1][j]) * mag_interactions[i][j][4]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);

                    b_field_z += (1*currents[0][j] + 0*currents[1][j]) * mag_interactions[i][j][5]*cexp(mag_interactions[i][j][7]*2*M_PI*freq*I/SPEED_OF_LIGHT);
                }
            }

            b_tot_squared = pow(creal(b_field_x),2) + pow(creal(b_field_y),2) + pow(creal(b_field_z),2) + pow(cimag(b_field_x),2) + pow(cimag(b_field_y),2) + pow(cimag(b_field_z),2);

            prox_current = M_PI*b_tot_squared*pow(freq*M_PI*conductivity*wire_radius*wire_radius,2);

            if((i == 0) || (i == num_points)){
                len = len_bases;
            }
            else{
                len = len_segs;
            }

            eq_resistance = len * (skin_loss_coeff + prox_current/mean_squared_coil_curr)/conductivity;
            fprintf(fp, "LD 4 0 %d 0 %lf 0\n", i+1, eq_resistance);
            //printf("\n%e", cabs(eq_resistance));
            if(i==num_points){
                //printf("\n%e", cabs(eq_resistance));
            }

        }
        (*num_iterations)++;


      }
      else{
        if((*num_iterations) >= max_iterations){
            printf("\nWARNING: Frequency step did not converge to the specified limits\nRel error: %e\n",rel_error);
        }
        else{
            printf("\nSuccessful frequency step in imposed limits using %d iterations\n",(*num_iterations));
        }
        return 1;
      }
    }
  }
  // Define a perfect ground plane (GN 1)
  fprintf(fp, "GN 1\n");

  // Use the extended kernel for wire
  fprintf(fp, "EK\n");

  // Excitation: Excite at the first segment (segment 1 of wire 1)
  fprintf(fp, "EX 0 1 1 0 1 0 \n");

  // Frequency definition (optional, using 1000 MHz as an example)
  fprintf(fp, "FR 0 1 0 0 %lf 0\n", freq/1e6);

  // Execute the simulation
  fprintf(fp, "XQ 0\n");

  // End of the NEC file
  fprintf(fp, "EN\n");

  fclose(fp);
  return 0;
}



// Function to find the byte position of the third line after the line where the given string is found
long find_string_in_file(const char *filename, const char *search_string, int lines_to_skip) {
  FILE *file;
  char *line = NULL;      // Pointer to dynamically allocated line buffer
  size_t len = 0;         // Holds the size of the buffer
  ssize_t read_len;       // Holds the length of the line read by getline()

  long position = 0;   // Position of the file pointer

  // Open the file in read mode
  file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return -1; // Return -1 if file can't be opened
  }

  // Read the file line by line using getline()
  while ((read_len = getline(&line, &len, file)) != -1) {

    // Check if the search string is in the current line
    if (strstr(line, search_string) != NULL) {
      //printf("line_found");
      // After finding the target string, skip the next two lines
      for (int i = 0; i < lines_to_skip; i++) {
        if ((read_len = getline(&line, &len, file)) == -1) {
          // If we reach EOF before skipping 2 lines, return -1
          free(line);
          fclose(file);
          return -1;
        }
      }

      // After skipping 2 lines, get the position of the third line
      position = ftell(file);
      free(line);  // Free the dynamically allocated buffer
      fclose(file);
      return position; // Return the byte position of the third line
    }
  }

  // Cleanup
  free(line);  // Free the dynamically allocated buffer
  fclose(file);

  // If the string is not found or file ends before reaching the third line, return -1
  return -1;
}



double complex read_impedance_data_nec_out(const char *filename, long offset) {
  // Open the file in read mode
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return -1;
  }

  // Move the file pointer to the specified offset
  if (fseek(file, offset, SEEK_SET) != 0) {
      perror("Error using fseek");
      fclose(file);
      return -1;
  }

  // Variables to hold the data
  int int1, int2;
  double d1, d2, d3, d4, d5, d6, d7, d8, d9;
  double complex reading =0;

  // Read the two integers and the rest as double values
  if (fscanf(file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&int1, &int2, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9) == 11) {
    reading = d5 + d6 * I;
  }
  else {
    printf("Error reading the data from the file at the specified offset.\n");
  }

  // Close the file
  fclose(file);

  return reading;
}

int read_simulation_currents_data(const char *filename, long offset, int num_segments, double complex *currents) {
  // Open the file in read mode
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return -1;
  }

  // Move the file pointer to the specified offset
  if (fseek(file, offset, SEEK_SET) != 0) {
      perror("Error using fseek");
      fclose(file);
      return -1;
  }

  // Variables to hold the data
  int i, int1, int2;
  double pos_x, pos_y, pos_z, len, curr_re, curr_im, curr_mag, curr_arg;

  for(i=0; i<num_segments; i++){

    // Read the two integers and the rest as double values
    if (fscanf(file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&int1, &int2, &pos_x, &pos_y, &pos_z, &len, &curr_re, &curr_im, &curr_mag, &curr_arg) == 10) {
      currents[i] = curr_re + curr_im * I;
    }
    else {
      printf("Error reading the data from the file at the specified offset.\n");
      return -1;
    }
  }
  // Close the file
  fclose(file);

  return 0;
}

// Function to append ".nec", ".out" and ".s1p" to the input string
void append_extensions(const char *input, char **output_nec, char **output_out, char **output_s_p) {
    // Calculate the length of the input string and the total length of the new strings
    size_t input_len = strlen(input);
    size_t nec_len = input_len + 4;  // for ".nec"
    size_t out_len = input_len + 4;  // for ".out"
    size_t s_p_len = input_len + 4;  // for ".s1p"

    // Allocate memory for the new strings (don't forget null terminators)
    *output_nec = (char *)malloc(nec_len + 1); // +1 for null terminator
    *output_out = (char *)malloc(out_len + 1); // +1 for null terminator
    *output_s_p = (char *)malloc(s_p_len + 1); // +1 for null terminator

    // Create the new strings by appending the extensions
    strcpy(*output_nec, input);   // Copy input to output_nec
    strcat(*output_nec, ".nec");  // Append ".nec"

    strcpy(*output_out, input);   // Copy input to output_out
    strcat(*output_out, ".out");  // Append ".out"

    strcpy(*output_s_p, input);   // Copy input to output_s_p
    strcat(*output_s_p, ".s1p");  // Append ".s1p"
}


void generate_touchstone(
    double complex *Z_data,   // Impedance data (complex array)
    double *freq_data,        // Frequency data (real array)
    int data_size,            // Size of the arrays
    char *filename,           // Output filename
    char mode,                // 'S' for S-parameters, 'Y' for Y-parameters, 'Z' for Z-parameters
    double Z0,                // Normalization impedance (Z0)
    char **comments           // List of comment lines (each line ends with \n), last pointer is NULL
) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error opening file %s\n", filename);
        return;
    }

    // Calculate the appropriate frequency unit based on geometric mean
    double scale;
    char *freq_unit = (char *)malloc(4 * sizeof(char));

    double freq_mean_exponent = log10(sqrt(freq_data[0] * freq_data[data_size - 1]));  // Geometric mean of the first and last frequencies

    if (freq_mean_exponent >= 7.5) {
        scale = 1e9;
        strcpy(freq_unit, "GHz");
    } else if (freq_mean_exponent >= 4.5) {
        scale = 1e6;
        strcpy(freq_unit, "MHz");
    } else if (freq_mean_exponent >= 1.5) {
        scale = 1e3;
        strcpy(freq_unit, "kHz");
    } else {
        scale = 1.0;
        strcpy(freq_unit, "Hz");
    }

    // Write comment lines to the file
    if (comments != NULL) {
        for (int i = 0; comments[i] != NULL; i++) {
            fprintf(file, "! %s", comments[i]);  // Write each comment line
        }
    }

    // Write the Touchstone file headers for version 2.0
    fprintf(file, "[Version] 2.0\n");
    fprintf(file, "# %s %c RI R %.2f\n", freq_unit, mode, Z0);  // Frequency units based on geometric mean scaling
    fprintf(file, "[Number of Ports] 1\n");
    fprintf(file, "[Number of Frequencies] %d\n", data_size);
    fprintf(file, "[Network Data]\n");

    // Convert the impedance (Z), admittance (Y), or scattering (S) parameters based on the mode
    for (int i = 0; i < data_size; i++) {
        double freq_scaled = freq_data[i] / scale;  // Scale frequency to the correct unit
        double complex Z = Z_data[i];              // Impedance at the current frequency

        if (mode == 'S') {
            // Calculate S11 from Z11
            double complex S11 = (Z - Z0) / (Z + Z0);
            fprintf(file, "%.6f %.6f %.6f\n", freq_scaled, creal(S11), cimag(S11));
        } else if (mode == 'Y') {
            // Calculate Y11 (admittance) from Z11
            double complex Y = 1.0 / Z;
            fprintf(file, "%.6f %.6f %.6f\n", freq_scaled, creal(Y), cimag(Y));
        } else if (mode == 'Z') {
            // Write the impedance directly
            fprintf(file, "%.6f %.6f %.6f\n", freq_scaled, creal(Z), cimag(Z));
        } else {
            printf("Unknown mode: %c\n", mode);
            fclose(file);
            return;
        }
    }

    // Close the file
    fclose(file);
    printf("Touchstone file created: %s\n", filename);
}

// Function to compute the single integral using Simpson's rule
double compute_prox_eff_geometry_integral(double R, double L, int n_intervals) {
    // n_intervals must be even for Simpson's rule
    if (n_intervals % 2 != 0) {
        n_intervals++;  // Increase by 1 if odd
    }

    double h = R / n_intervals;
    double integral = 0.0;

    for (int i = 0; i <= n_intervals; i++) {
        double x = i * h;
        // Ensure the square root argument is non-negative (it should be by definition)
        double sqrt_term = sqrt(R*R - x*x);
        double f = (x*x * sqrt_term) / (2*x + L);

        if (i == 0 || i == n_intervals) {
            integral += f;
        } else if (i % 2 == 1) {
            integral += 4 * f;
        } else {
            integral += 2 * f;
        }
    }

    integral *= h / 3.0;
    return integral;
}

int main() {
    // Declare variables for user input
    Point3D* point_array = NULL;
    double wire_rad, coil_rad, coil_length, freq, height, freq_step, conductivity = 59.6e6, segment_length;
    int n_turns, freq_step_num, n_points, i=0, iterations = 0, max_iterations = 300;
    char file_name[100];
    const char* impedance_line_preamble = "                        --------- ANTENNA INPUT PARAMETERS ---------";
    const char* currents_line_preamble = "                           -------- CURRENTS AND LOCATION --------";
    char* input_filename = NULL;
    char* out_filename = NULL;
    char* prompt_base = "nec2c -i ";
    char* prompt_mid = " -o ";
    char* full_prompt = NULL;
    long impedance_line_offset;
    long currents_line_offset;
    //double complex result = 0;
    char** comment_lines= NULL;
    double complex* impedance_array = NULL;
    double complex* currents[2];
    double* frequency_array = NULL;
    char* S_parameter_filename = NULL;

    int use_skin_effect = 1;
    int use_proximity_effect = 1;


    // Prompt the user for input values
    printf("Enter the wire radius (in meters, e.g., 1e-4): ");
    scanf("%lf", &wire_rad);

    printf("Enter the coil radius (in meters, e.g., 2.5e-3): ");
    scanf("%lf", &coil_rad);

    printf("Enter the coil length (in meters, e.g., 1e-2): ");
    scanf("%lf", &coil_length);

    printf("Enter the number of turns (e.g., 4): ");
    scanf("%d", &n_turns);

    printf("Enter the starting frequency (in Hz, e.g., 300e6): ");
    scanf("%lf", &freq);

    printf("Enter the extra height (in meters, e.g., 1.6e-3): ");
    scanf("%lf", &height);

    printf("Enter the frequency step (in Hz, e.g., 100e6): ");
    scanf("%lf", &freq_step);

    printf("Enter the number of frequency steps (e.g., 58): ");
    scanf("%d", &freq_step_num);

    printf("Enter the output NEC file name (e.g., first_coil): ");
    scanf("%s", file_name);

    if(freq_step_num <1){
        return 0;
    }

    // Reserve memory for simulation results
    impedance_array = (double complex*)malloc(freq_step_num * sizeof(double complex));
    frequency_array = (double*)malloc(freq_step_num * sizeof(double));

    // Create the filenames for input and output files
    append_extensions(file_name, &input_filename, &out_filename, &S_parameter_filename);

    // Generate the command for calling nec2c with the coil data
    size_t command_length = strlen(prompt_base) + strlen(input_filename) + strlen(prompt_mid)+ strlen(out_filename);
    full_prompt = (char *)malloc(command_length + 1);
    strcpy(full_prompt, prompt_base);
    strcat(full_prompt, input_filename);
    strcat(full_prompt, prompt_mid);
    strcat(full_prompt, out_filename);

    // Generate coil points based on user input
    point_array = generate_coil_points(coil_rad, coil_length, n_turns, wire_rad, freq, freq + freq_step_num * freq_step, height, &n_points, &segment_length, 0, 16);

    double proximity_effect_constant, proximity_effect_constant_bases;
    double (*(*magnetic_interactions))[8] = NULL;
    if(use_proximity_effect){

        proximity_effect_constant = -4 * compute_prox_eff_geometry_integral(wire_rad, segment_length, 1000) * pow(segment_length * conductivity, 2);
        proximity_effect_constant_bases = -4 * compute_prox_eff_geometry_integral(wire_rad, height, 1000) * pow(height * conductivity, 2);

        printf("\n%e %e %e\n",compute_prox_eff_geometry_integral(wire_rad, segment_length, 1000), segment_length, proximity_effect_constant);

        double dist_x, dist_y, dist_z, dist_z_ref, dl_x, dl_y, dl_z, dl_z_ref, dot_prod, dot_prod_ref ;

        magnetic_interactions = (double (*(*))[8])malloc(sizeof(double(*)[8])*(n_points + 1));
        magnetic_interactions[0] = (double (*)[8])malloc(sizeof(double[8])*(n_points + 1));

        for(i = 0; i<7; i++){
            magnetic_interactions[0][0][i] = 0;
        }
        magnetic_interactions[0][0][7] = height;

        for(i=1; i<n_points; i++){
            dist_x= 0.5*(point_array[i-1].x + point_array[i].x) - point_array[0].x;
            dist_y= 0.5*(point_array[i-1].y + point_array[i].y) - point_array[0].y;
            dist_z= 0.5*(point_array[i-1].z + point_array[i].z) - 0.5*height;
            dist_z_ref = -0.5*(point_array[i-1].z + point_array[i].z) - 0.5*height;

            dl_x = point_array[i].x - point_array[i-1].x;
            dl_y = point_array[i].y - point_array[i-1].y;
            dl_z = point_array[i].z - point_array[i-1].z;
            dl_z_ref = dl_z * -1;

            magnetic_interactions[0][i][6] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
            magnetic_interactions[0][i][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

            magnetic_interactions[0][i][0] = dist_y*dl_z - dl_y*dist_z;
            magnetic_interactions[0][i][1]= dist_z*dl_x - dl_z*dist_x;
            magnetic_interactions[0][i][2] = dist_x*dl_y - dl_x*dist_y;

            dot_prod = point_array[0].z*magnetic_interactions[0][i][2];

            magnetic_interactions[0][i][0] = MU_0*(magnetic_interactions[0][i][0] - (dot_prod*point_array[0].x)/pow(point_array[0].z,2))/(pow(magnetic_interactions[0][i][6],3)*4*M_PI);
            magnetic_interactions[0][i][1] = MU_0*(magnetic_interactions[0][i][1] - (dot_prod*point_array[0].y)/pow(point_array[0].z,2))/(pow(magnetic_interactions[0][i][6],3)*4*M_PI);
            magnetic_interactions[0][i][2] = 0;

            magnetic_interactions[0][i][3] = dist_y*dl_z_ref - dl_y*dist_z_ref;
            magnetic_interactions[0][i][4] = dist_z_ref*dl_x - dl_z_ref*dist_x;
            magnetic_interactions[0][i][5] = 0;

            dot_prod_ref = dot_prod;

            magnetic_interactions[0][i][3] = MU_0*(magnetic_interactions[0][i][3] - (dot_prod_ref*point_array[0].x)/pow(point_array[0].z,2))/(pow(magnetic_interactions[0][i][7],3)*4*M_PI);
            magnetic_interactions[0][i][4] = MU_0*(magnetic_interactions[0][i][4] - (dot_prod_ref*point_array[0].y)/pow(point_array[0].z,2))/(pow(magnetic_interactions[0][i][7],3)*4*M_PI);
            magnetic_interactions[0][i][5] = 0;
        }

        //añadir interaccion de [n_points] (shorted element)
        dist_x = point_array[n_points-1].x - point_array[0].x;
        dist_y = point_array[n_points-1].y - point_array[0].y;
        dist_z = 0;
        dist_z_ref = -1 * height;

        dl_x = 0;
        dl_y = 0;
        dl_z = -1 * point_array[n_points-1].z;
        dl_z_ref = dl_z * -1;

        magnetic_interactions[0][n_points][6] = sqrt(pow(dist_x,2) + pow(dist_y,2));
        magnetic_interactions[0][n_points][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

        magnetic_interactions[0][n_points][0] = MU_0*dist_y*dl_z/(pow(magnetic_interactions[0][n_points][6],3)*4*M_PI);
        magnetic_interactions[0][n_points][1]= -MU_0*dl_z*dist_x/(pow(magnetic_interactions[0][n_points][6],3)*4*M_PI);
        magnetic_interactions[0][n_points][2] = 0;

        dot_prod = 0;

        magnetic_interactions[0][n_points][3] = MU_0*dist_y*dl_z_ref/(pow(magnetic_interactions[0][n_points][7],3)*4*M_PI);
        magnetic_interactions[0][n_points][4] = -MU_0*dl_z_ref*dist_x/(pow(magnetic_interactions[0][n_points][7],3)*4*M_PI);
        magnetic_interactions[0][n_points][5] = 0;

        dot_prod_ref = 0;

        for(i = 1; i<(n_points); i++){
            magnetic_interactions[i] = (double(*)[8])malloc(sizeof(double[8])*(n_points + 1));

            dist_x = point_array[0].x -0.5*(point_array[i-1].x + point_array[i].x);
            dist_y = point_array[0].y -0.5*(point_array[i-1].x + point_array[i].x);
            dist_z = 0.5*(height - (point_array[i-1].z + point_array[i].z));
            dist_z_ref = -0.5*(height + point_array[i-1].z + point_array[i].z);

            dl_x = point_array[0].x;
            dl_y = point_array[0].y;
            dl_z = point_array[0].z;
            dl_z_ref = dl_z * -1;

            magnetic_interactions[i][0][6] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
            magnetic_interactions[i][0][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

            magnetic_interactions[i][0][0] = dist_y*dl_z - dl_y*dist_z;
            magnetic_interactions[i][0][1]= dist_z*dl_x - dl_z*dist_x;
            magnetic_interactions[i][0][2] = dist_x*dl_y - dl_x*dist_y;

            dot_prod = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][0][0] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][0][1] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][0][2];

            magnetic_interactions[i][0][0] = MU_0*(magnetic_interactions[i][0][0] - (dot_prod*point_array[0].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][6],3)*4*M_PI);
            magnetic_interactions[i][0][1] = MU_0*(magnetic_interactions[i][0][1] - (dot_prod*point_array[0].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][6],3)*4*M_PI);
            magnetic_interactions[i][0][2] = MU_0*(magnetic_interactions[i][0][2] - (dot_prod*point_array[0].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][6],3)*4*M_PI);

            magnetic_interactions[i][0][3] = dist_y*dl_z_ref - dl_y*dist_z_ref;
            magnetic_interactions[i][0][4] = dist_z_ref*dl_x - dl_z_ref*dist_x;
            magnetic_interactions[i][0][5] = dist_x*dl_y - dl_x*dist_y;

            dot_prod_ref = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][0][3] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][0][4] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][0][5];

            magnetic_interactions[i][0][3] = MU_0*(magnetic_interactions[i][0][3] - dot_prod_ref*(point_array[i].x - point_array[i-1].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][7],3)*4*M_PI);
            magnetic_interactions[i][0][4] = MU_0*(magnetic_interactions[i][0][4] - dot_prod_ref*(point_array[i].y - point_array[i-1].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][7],3)*4*M_PI);
            magnetic_interactions[i][0][5] = MU_0*(magnetic_interactions[i][0][5] - dot_prod_ref*(point_array[i].z - point_array[i-1].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][0][7],3)*4*M_PI);

            for(int j = 1; j<n_points; j++){
                if(i != j){
                    //arreglar para caso general
                    dist_x = 0.5*((point_array[j-1].x + point_array[j].x) - (point_array[i-1].x + point_array[i].x));
                    dist_y = 0.5*((point_array[j-1].y + point_array[j].y) - (point_array[i-1].y + point_array[i].y));
                    dist_z = 0.5*((point_array[j-1].z + point_array[j].z) - (point_array[i-1].z + point_array[i].z));
                    dist_z_ref = -0.5*(point_array[j-1].z + point_array[j].z + point_array[i-1].z + point_array[i].z);

                    dl_x = point_array[j].x - point_array[j-1].x;
                    dl_y = point_array[j].y - point_array[j-1].y;
                    dl_z = point_array[j].z - point_array[j-1].z;
                    dl_z_ref = dl_z * -1;

                    magnetic_interactions[i][j][6] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
                    magnetic_interactions[i][j][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

                    magnetic_interactions[i][j][0] = dist_y*dl_z - dl_y*dist_z;
                    magnetic_interactions[i][j][1]= dist_z*dl_x - dl_z*dist_x;
                    magnetic_interactions[i][j][2] = dist_x*dl_y - dl_x*dist_y;

                    dot_prod = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][j][0] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][j][1] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][j][2];

                    magnetic_interactions[i][j][0] = MU_0*(magnetic_interactions[i][j][0] - dot_prod*(point_array[i].x - point_array[i-1].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][6],3)*4*M_PI);
                    magnetic_interactions[i][j][1] = MU_0*(magnetic_interactions[i][j][1] - dot_prod*(point_array[i].y - point_array[i-1].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][6],3)*4*M_PI);
                    magnetic_interactions[i][j][2] = MU_0*(magnetic_interactions[i][j][2] - dot_prod*(point_array[i].z - point_array[i-1].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][6],3)*4*M_PI);

                    magnetic_interactions[i][j][3] = dist_y*dl_z_ref - dl_y*dist_z_ref;
                    magnetic_interactions[i][j][4] = dist_z_ref*dl_x - dl_z_ref*dist_x;
                    magnetic_interactions[i][j][5] = dist_x*dl_y - dl_x*dist_y;

                    dot_prod_ref = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][j][3] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][j][4] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][j][5];

                    magnetic_interactions[i][j][3] = MU_0*(magnetic_interactions[i][j][3] - dot_prod_ref*(point_array[i].x - point_array[i-1].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][7],3)*4*M_PI);
                    magnetic_interactions[i][j][4] = MU_0*(magnetic_interactions[i][j][4] - dot_prod_ref*(point_array[i].y - point_array[i-1].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][7],3)*4*M_PI);
                    magnetic_interactions[i][j][5] = MU_0*(magnetic_interactions[i][j][5] - dot_prod_ref*(point_array[i].z - point_array[i-1].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][7],3)*4*M_PI);
                }
                else{
                    //Hacer caso de interaccion de propio elemento
                    dist_x = 0;
                    dist_y = 0;
                    dist_z_ref = -(point_array[j-1].z + point_array[j].z);

                    dl_x = point_array[j].x - point_array[j-1].x;
                    dl_y = point_array[j].y - point_array[j-1].y;
                    dl_z_ref = point_array[j-1].z - point_array[j].z;

                    magnetic_interactions[i][j][6] = 0;
                    magnetic_interactions[i][j][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

                    magnetic_interactions[i][j][0] = 0;
                    magnetic_interactions[i][j][1] = 0;
                    magnetic_interactions[i][j][2] = 0;

                    magnetic_interactions[i][j][3] = -1 * dl_y*dist_z_ref;
                    magnetic_interactions[i][j][4] = dist_z_ref*dl_x;
                    magnetic_interactions[i][j][5] = 0;

                    dot_prod_ref = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][j][3] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][j][4];

                    magnetic_interactions[i][j][3] = MU_0*(magnetic_interactions[i][j][3] - dot_prod_ref*(point_array[i].x - point_array[i-1].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][7],3)*4*M_PI);
                    magnetic_interactions[i][j][4] = MU_0*(magnetic_interactions[i][j][4] - dot_prod_ref*(point_array[i].y - point_array[i-1].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][j][7],3)*4*M_PI);
                    magnetic_interactions[i][j][5] = 0;
                }

            }
            //elemnto final [n_points]
            dist_x = point_array[n_points-1].x -0.5*(point_array[i-1].x + point_array[i].x);
            dist_y = point_array[n_points-1].y -0.5*(point_array[i-1].x + point_array[i].x);
            dist_z = 0.5*(height - (point_array[i-1].z + point_array[i].z));
            dist_z_ref = -0.5*(height + point_array[i-1].z + point_array[i].z);

            dl_x = point_array[n_points-1].x;
            dl_y = point_array[n_points-1].y;
            dl_z = point_array[n_points-1].z;
            dl_z_ref = dl_z * -1;

            magnetic_interactions[i][n_points][6] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
            magnetic_interactions[i][n_points][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

            magnetic_interactions[i][n_points][0] = dist_y*dl_z - dl_y*dist_z;
            magnetic_interactions[i][n_points][1]= dist_z*dl_x - dl_z*dist_x;
            magnetic_interactions[i][n_points][2] = dist_x*dl_y - dl_x*dist_y;

            dot_prod = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][n_points][0] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][n_points][1] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][n_points][2];

            magnetic_interactions[i][n_points][0] = MU_0*(magnetic_interactions[i][n_points][0] - (dot_prod*point_array[0].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][6],3)*4*M_PI);
            magnetic_interactions[i][n_points][1] = MU_0*(magnetic_interactions[i][n_points][1] - (dot_prod*point_array[0].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][6],3)*4*M_PI);
            magnetic_interactions[i][n_points][2] = MU_0*(magnetic_interactions[i][n_points][2] - (dot_prod*point_array[0].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][6],3)*4*M_PI);

            magnetic_interactions[i][n_points][3] = dist_y*dl_z_ref - dl_y*dist_z_ref;
            magnetic_interactions[i][n_points][4] = dist_z_ref*dl_x - dl_z_ref*dist_x;
            magnetic_interactions[i][n_points][5] = dist_x*dl_y - dl_x*dist_y;

            dot_prod_ref = (point_array[i].x - point_array[i-1].x)*magnetic_interactions[i][n_points][3] + (point_array[i].y - point_array[i-1].y)*magnetic_interactions[i][n_points][4] + (point_array[i].z - point_array[i-1].z)*magnetic_interactions[i][n_points][5];

            magnetic_interactions[i][n_points][3] = MU_0*(magnetic_interactions[i][n_points][3] - dot_prod_ref*(point_array[i].x - point_array[i-1].x)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][7],3)*4*M_PI);
            magnetic_interactions[i][n_points][4] = MU_0*(magnetic_interactions[i][n_points][4] - dot_prod_ref*(point_array[i].y - point_array[i-1].y)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][7],3)*4*M_PI);
            magnetic_interactions[i][n_points][5] = MU_0*(magnetic_interactions[i][n_points][5] - dot_prod_ref*(point_array[i].z - point_array[i-1].z)/pow(segment_length,2))/(pow(magnetic_interactions[i][n_points][7],3)*4*M_PI);
        }
        magnetic_interactions[n_points] = (double (*)[8])malloc(sizeof(double[8])*(n_points + 1));

        //añadir interaccion de [n_points] (shorted element)
        dist_x = point_array[0].x - point_array[n_points-1].x;
        dist_y = point_array[0].y - point_array[n_points-1].y;
        dist_z = 0;
        dist_z_ref = -1 * height;

        dl_x = 0;
        dl_y = 0;
        dl_z = -1 * point_array[n_points-1].z;
        dl_z_ref = dl_z * -1;

        magnetic_interactions[n_points][0][6] = sqrt(pow(dist_x,2) + pow(dist_y,2));
        magnetic_interactions[n_points][0][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

        magnetic_interactions[n_points][0][0] = MU_0*dist_y*dl_z/(pow(magnetic_interactions[n_points][0][6],3)*4*M_PI);
        magnetic_interactions[n_points][0][1]= -MU_0*dl_z*dist_x/(pow(magnetic_interactions[n_points][0][6],3)*4*M_PI);
        magnetic_interactions[n_points][0][2] = 0;

        dot_prod = 0;

        magnetic_interactions[n_points][0][3] = MU_0*dist_y*dl_z_ref/(pow(magnetic_interactions[n_points][0][7],3)*4*M_PI);
        magnetic_interactions[n_points][0][4] = -MU_0*dl_z_ref*dist_x/(pow(magnetic_interactions[n_points][0][7],3)*4*M_PI);
        magnetic_interactions[n_points][0][5] = 0;

        dot_prod_ref = 0;

        for(i=1; i<n_points;i++){
            dist_x= 0.5*(point_array[i-1].x + point_array[i].x) - point_array[n_points-1].x;
            dist_y= 0.5*(point_array[i-1].y + point_array[i].y) - point_array[n_points-1].y;
            dist_z= 0.5*(point_array[i-1].z + point_array[i].z) - 0.5*height;
            dist_z_ref = -0.5*(point_array[i-1].z + point_array[i].z) - 0.5*height;

            dl_x = point_array[i].x - point_array[i-1].x;
            dl_y = point_array[i].y - point_array[i-1].y;
            dl_z = point_array[i].z - point_array[i-1].z;
            dl_z_ref = dl_z * -1;

            magnetic_interactions[n_points][i][6] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
            magnetic_interactions[n_points][i][7] = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z_ref,2));

            magnetic_interactions[n_points][i][0] = dist_y*dl_z - dl_y*dist_z;
            magnetic_interactions[n_points][i][1]= dist_z*dl_x - dl_z*dist_x;
            magnetic_interactions[n_points][i][2] = dist_x*dl_y - dl_x*dist_y;

            dot_prod = point_array[n_points-1].z*magnetic_interactions[n_points][i][2];

            magnetic_interactions[n_points][i][0] = MU_0*(magnetic_interactions[n_points][i][0] - (dot_prod*point_array[n_points-1].x)/pow(point_array[n_points-1].z,2))/(pow(magnetic_interactions[n_points][i][6],3)*4*M_PI);
            magnetic_interactions[n_points][i][1] = MU_0*(magnetic_interactions[n_points][i][1] - (dot_prod*point_array[n_points-1].y)/pow(point_array[n_points-1].z,2))/(pow(magnetic_interactions[n_points][i][6],3)*4*M_PI);
            magnetic_interactions[n_points][i][2] = 0;

            magnetic_interactions[n_points][i][3] = dist_y*dl_z_ref - dl_y*dist_z_ref;
            magnetic_interactions[n_points][i][4] = dist_z_ref*dl_x - dl_z_ref*dist_x;
            magnetic_interactions[n_points][i][5] = 0;

            dot_prod_ref = dot_prod;

            magnetic_interactions[n_points][i][3] = MU_0*(magnetic_interactions[n_points][i][3] - (dot_prod_ref*point_array[n_points-1].x)/pow(point_array[n_points-1].z,2))/(pow(magnetic_interactions[n_points][i][7],3)*4*M_PI);
            magnetic_interactions[n_points][i][4] = MU_0*(magnetic_interactions[n_points][i][4] - (dot_prod_ref*point_array[n_points-1].y)/pow(point_array[n_points-1].z,2))/(pow(magnetic_interactions[n_points][i][7],3)*4*M_PI);
            magnetic_interactions[n_points][i][5] = 0;
        }

        for(i=0; i<7; i++){
            magnetic_interactions[n_points][n_points][i] = 0;
        }
        magnetic_interactions[n_points][n_points][7] = height;
    }


    for(int i = 0; i <= n_points; i++) {
        for(int j = 0; j <= n_points; j++) {
            //printf("Element [%d][%d]: ", i, j);
            for(int k = 0; k < 8; k++) {
                //printf("%lf ", magnetic_interactions[i][j][k]);
            }
            //printf("\n");
        }
        //printf("\n");
    }
    //printf("String content: %s\n", out_filename);
    if(use_proximity_effect){
        currents[0] = (double complex*)malloc((n_points+1) * sizeof(double complex));
        currents[1] = (double complex*)malloc((n_points+1) * sizeof(double complex));
        if(currents[0] == NULL || currents[1] == NULL){
            printf("\nMemory allocation error with currents array\n");
            return 1;
        }


        iterations = 1;

        //Create a nec2 simulation with only the skin effect only resistance in the wires
        create_nec_file(input_filename, point_array, n_points, wire_rad, freq, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,2,&iterations,max_iterations,0.001*(n_points+1));

        // Execute the simulation using nec2c engine
        system(full_prompt);
        //printf("\nok\n");

        //Extract the currents in the simulation and update the currents array
        currents_line_offset = find_string_in_file(out_filename, currents_line_preamble, 4);
        read_simulation_currents_data(out_filename,currents_line_offset,n_points+1,currents[0]);

        // Create the NEC file with the user inputs and iterate to find the solution
        while(1 != create_nec_file(input_filename, point_array, n_points, wire_rad, freq, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,3,&iterations,max_iterations,0.001*(n_points+1))){
            // Execute the simulation using nec2c engine
            system(full_prompt);
            //printf("\nok\n");
            //Extract the currents in the simulation and update the currents array

            for(int i = 0; i <= n_points; i++){
                currents[1][i]= currents[0][i];
            }

            currents_line_offset = find_string_in_file(out_filename, currents_line_preamble, 4);

            read_simulation_currents_data(out_filename,currents_line_offset,n_points+1,currents[0]);
        }

        // Measure impedance parameter from output file
        impedance_line_offset = find_string_in_file(out_filename, impedance_line_preamble, 2);

        impedance_array[0] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
        frequency_array[0] = freq;

        for(i = 1; i<freq_step_num; i++){
            iterations = 1;

            //Create a nec2 simulation with only the skin effect only resistance in the wires
            create_nec_file(input_filename, point_array, n_points, wire_rad, freq + i*freq_step, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,2,&iterations,max_iterations,0.001*(n_points+1));

            // Execute the simulation using nec2c engine
            system(full_prompt);
            //printf("\nok\n");

            //Extract the currents in the simulation and update the currents array
            currents_line_offset = find_string_in_file(out_filename, currents_line_preamble, 4);
            read_simulation_currents_data(out_filename,currents_line_offset,n_points+1,currents[0]);

            while(1 != create_nec_file(input_filename, point_array, n_points, wire_rad, freq + i*freq_step, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,3,&iterations,max_iterations,0.001*(n_points+1))){
                // Execute the simulation using nec2c engine
                system(full_prompt);
                //printf("\nok\n");
                //Extract the currents in the simulation and update the currents array
                for(int i = 0; i <= n_points; i++){
                    currents[1][i]= currents[0][i];
                }
                currents_line_offset = find_string_in_file(out_filename, currents_line_preamble, 4);
                read_simulation_currents_data(out_filename,currents_line_offset,n_points+1,currents[0]);

            }

            // Measure impedance parameter from output file
            impedance_array[i] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
            frequency_array[i] = freq + i*freq_step;
            //printf("Freq(%lf) Re(%lf) Im(%lf)\n", freq +  i*freq_step, creal(impedance_array[i]), cimag(impedance_array[i]));
        }
    }
    else{
        if(use_skin_effect){
            // Create the NEC file with the user inputs
            create_nec_file(input_filename, point_array, n_points, wire_rad, freq, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,2,&iterations,max_iterations,0.01/(n_points+1));

            // Execute the simulation using nec2c engine
            system(full_prompt);

            // Measure impedance parameter from output file
            impedance_line_offset = find_string_in_file(out_filename, impedance_line_preamble, 2);

            impedance_array[0] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
            frequency_array[0] = freq;
            for(i = 1; i<freq_step_num; i++){

                // Create the NEC file with the user inputs
                create_nec_file(input_filename, point_array, n_points, wire_rad, freq + i*freq_step, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,2,&iterations,max_iterations,0.01/(n_points+1));

                // Execute the simulation using nec2c engine
                system(full_prompt);

                // Measure impedance parameter from output file
                impedance_array[i] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
                frequency_array[i] = freq + i*freq_step;
                //printf("Freq(%lf) Re(%lf) Im(%lf)\n", freq +  i*freq_step, creal(impedance_array[i]), cimag(impedance_array[i]));
            }
        }
        else{
            // Create the NEC file with the user inputs
            create_nec_file(input_filename, point_array, n_points, wire_rad, freq, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,1,&iterations,max_iterations,0.01/(n_points+1));

            // Execute the simulation using nec2c engine
            system(full_prompt);

            // Measure impedance parameter from output file
            impedance_line_offset = find_string_in_file(out_filename, impedance_line_preamble, 2);

            impedance_array[0] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
            frequency_array[0] = freq;

            for(i = 1; i<freq_step_num; i++){

                // Create the NEC file with the user inputs
                create_nec_file(input_filename, point_array, n_points, wire_rad, freq + i*freq_step, conductivity,currents,magnetic_interactions,proximity_effect_constant, proximity_effect_constant_bases,1,&iterations,max_iterations,0.01/(n_points+1));

                // Execute the simulation using nec2c engine
                system(full_prompt);

                // Measure impedance parameter from output file
                impedance_array[i] = read_impedance_data_nec_out(out_filename, impedance_line_offset);
                frequency_array[i] = freq + i*freq_step;
                //printf("Freq(%lf) Re(%lf) Im(%lf)\n", freq +  i*freq_step, creal(impedance_array[i]), cimag(impedance_array[i]));
            }
        }
    }
    printf("\n%lf %lf %lf\n", point_array[n_points-1].x, point_array[n_points-1].y, point_array[n_points-1].z);
    printf("\nFINISHED SIMULATING\n");

    generate_touchstone(
        impedance_array,   // Impedance data (complex array)
        frequency_array,        // Frequency data (real array)
        freq_step_num,            // Size of the arrays
        S_parameter_filename,           // Output filename
        'S',                // 'S' for S-parameters, 'Y' for Y-parameters, 'Z' for Z-parameters
        50,                // Normalization impedance (Z0)
        comment_lines          // List of comment lines (each line ends with \n), last pointer is NULL
    );

    // Free allocated memory for point_array
    free(point_array);
    free(input_filename);
    free(out_filename);
    free(full_prompt);
    free(impedance_array);
    free(frequency_array);

    return 0;
}

