

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#define FAILURE 0
#define SUCCESS !FAILURE

#define USER_NAME "acr18by"		//replace with your user name

typedef enum MODE { CPU, OPENMP, CUDA, ALL } MODE;
typedef enum OUTPUT_MODE { PPM_BINARY, PPM_PLAIN_TEXT } OUTPUT_MODE;

void print_help();
int process_command_line(int argc, char *argv[]);
unsigned char ** read_data(const char* fname);
unsigned char * gpu_cal(unsigned char *gpu_data);
unsigned char * gpu_cal_optimised(unsigned char *gpu_data);
void cpu_cal();
void openmp_cal();
int output(char * fname);

int c = 0;
unsigned int width = 0;
unsigned int height = 0;
unsigned char ** data;
char *in_file;
char ftype[2];
char *out_file;
int r, g, b;

MODE execution_mode = CPU;
OUTPUT_MODE output_mode = PPM_BINARY;


int main(int argc, char *argv[])
{
	if (process_command_line(argc, argv) == FAILURE)
		return 1;


	//TODO: read input image file (either binary or plain text PPM) 
	printf("Reading data from %s \n", in_file);
	data = read_data(in_file);

	//TODO: execute the mosaic filter based on the mode
	switch (execution_mode) {
	case (CPU): {
		// TODO: starting timing here
		clock_t start = clock(), diff;
		// TODO: calculate the average colour value

		cpu_cal();

		// Output the average colour value for the image
		printf("CPU Average image colour red = %d, green = %d, blue = %d \n", r, g, b);
		// TODO: end timing here
		diff = clock() - start;
		int msec = diff * 1000 / CLOCKS_PER_SEC;
		printf("CPU mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);


		break;
	}
	case (OPENMP): {
		//TODO: starting timing here
		//clock_t start = clock(), diff;
		//TODO: calculate the average colour value
		//double begin, diff;
		//begin = omp_get_wtime();
		openmp_cal();
		// Output the average colour value for the image
		printf("OPENMP Average image colour red = %d, green = %d, blue = %d \n", r, g, b);

		////TODO: end timing here
		//diff = omp_get_wtime() - begin;
		//int msec = diff * 1000;
		//printf("OPENMP mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);
		break;
	}
	case (CUDA): {
		cudaEvent_t start, stop;
		float milliseconds = 0;

		unsigned char *gpu_data;
		size_t size = height * width * 3 * sizeof(unsigned char);
		gpu_data = (unsigned char *)malloc(size);
		// transfer data from 2d to 1d
		for (int i = 0, i_1d = 0; i < height; i++) {
			for (int j = 0; j < width * 3; j++, i_1d++) {
				*(gpu_data + i_1d) = *(*(data + i) + j);
			}
		}

		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start);
		gpu_data = gpu_cal(gpu_data);
		//gpu_data = gpu_cal_optimised(gpu_data);
		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&milliseconds, start, stop);
		printf("CUDA Average image colour red = %d, green = %d, blue = %d \n", r, g, b);
		printf("Execution time is %f ms\n", milliseconds);

		//transfer data from 1d to 2d for output
		for (int i = 0, i_1d = 0; i < height; i++) {
			for (int j = 0; j < width * 3; j++, i_1d++) {
				*(*(data + i) + j) = *(gpu_data + i_1d);
			}
		}
		free(gpu_data);

		output(out_file);
		break;
	}
	case (ALL): {
		//TODO
		clock_t start = clock(), diff;
		cudaEvent_t c_start, c_stop;
		float milliseconds = 0;


		// CPU MODE
		cpu_cal();
		printf("\nCPU Average image colour red = %d, green = %d, blue = %d \n", r, g, b);
		diff = clock() - start;
		int msec = diff * 1000 / CLOCKS_PER_SEC;
		printf("CPU mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);
		start = clock();
		
		// OPENMP MODE
		openmp_cal();
		printf("\nOPENMP Average image colour red = %d, green = %d, blue = %d \n", r, g, b);
		diff = clock() - start;
		msec = diff * 1000 / CLOCKS_PER_SEC;
		printf("OPENMP mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);

		// CUDA MODE
		unsigned char *gpu_data;
		size_t size = height * width * 3 * sizeof(unsigned char);
		gpu_data = (unsigned char *)malloc(size);
		// transfer data from 2d to 1d
		for (int i = 0, i_1d = 0; i < height; i++) {
			for (int j = 0; j < width * 3; j++, i_1d++) {
				*(gpu_data + i_1d) = *(*(data + i) + j);
			}
		}
		// CUDA TIME START HERE
		cudaEventCreate(&c_start);
		cudaEventCreate(&c_stop);
		cudaEventRecord(c_start);

		gpu_data = gpu_cal(gpu_data);
		//gpu_data = gpu_cal_optimised(gpu_data);
		cudaEventRecord(c_stop);
		cudaEventSynchronize(c_stop);
		cudaEventElapsedTime(&milliseconds, c_start, c_stop);
		printf("\nCUDA mode execution time took %d s and %dms\n", (int)milliseconds / 1000, (int)milliseconds % 1000);
		printf("CUDA Average image colour red = %d, green = %d, blue = %d \n", r, g, b);

		//transfer data from 1d to 2d for output
		for (int i = 0, i_1d = 0; i < height; i++) {
			for (int j = 0; j < width * 3; j++, i_1d++) {
				*(*(data + i) + j) = *(gpu_data + i_1d);
			}
		}
		free(gpu_data);

		output(out_file);


		break;
	}
	}

	free(data);
	getchar();

	return 0;
}

int process_command_line(int argc, char *argv[]) {
	if (argc < 7) {
		fprintf(stderr, "Error: Missing program arguments. Correct usage is...\n");
		print_help();
		return FAILURE;
	}
	//first argument is always the executable name

	//read in the non optional command line arguments

	c = atoi(argv[1]);

	if (c <= 0) {
		printf("The value of c is invalid.");
		return FAILURE;

	}

	c = pow(2.0, (double)(int)log2(c)); // change the value of c to be valid


	if (!strcmp(argv[2], "CPU")) { execution_mode = CPU; };
	if (!strcmp(argv[2], "OPENMP")) { execution_mode = OPENMP; };
	if (!strcmp(argv[2], "CUDA")) { execution_mode = CUDA; };
	if (!strcmp(argv[2], "ALL")) { execution_mode = ALL; };
	//TODO: read in the input image name
	in_file = argv[4];
	//TODO: read in the output image name
	out_file = argv[6];
	//TODO: read in any optional part 3 arguments
	if (argc > 8) {
		if (!strcmp(argv[8], "PPM_BINARY")) { output_mode = PPM_BINARY; };
		if (!strcmp(argv[8], "PPM_PLAIN_TEXT")) { output_mode = PPM_PLAIN_TEXT; };
	}
			

	

	return SUCCESS;
}

void print_help() {
	printf("mosaic_%s C M -i input_file -o output_file [options]\n", USER_NAME);

	printf("where:\n");
	printf("\tC              Is the mosaic cell size which should be any positive\n"
		"\t               power of 2 number \n");
	printf("\tM              Is the mode with a value of either CPU, OPENMP, CUDA or\n"
		"\t               ALL. The mode specifies which version of the simulation\n"
		"\t               code should execute. ALL should execute each mode in\n"
		"\t               turn.\n");
	printf("\t-i input_file  Specifies an input image file\n");
	printf("\t-o output_file Specifies an output image file which will be used\n"
		"\t               to write the mosaic image\n");
	printf("[options]:\n");
	printf("\t-f ppm_format  PPM image output format either PPM_BINARY (default) or \n"
		"\t               PPM_PLAIN_TEXT\n ");
}

/* Read header information of the file*/
FILE *read_header(FILE *fp) {
	char read_line[10];
	while (1) {
		// exit if reading to the end of file
		if (fgets(read_line, sizeof(read_line), fp) == NULL) {
			return FAILURE;
		}
		// exit if reading to the end line of header
		if (strncmp(read_line, "255", 3) == 0) {
			//size = str_cat(size, input);
			break;
		}
		// file format (either P3 or P6)
		if (strncmp(read_line, "P3", 2) == 0) {
			strcpy(ftype, "P3");
		}
		else if (strncmp(read_line, "P6", 2) == 0) {
			strcpy(ftype, "P6");
		}
		// skip if reading to command line
		else if (strncmp(read_line, "#", 1) == 0) {
			continue;
		}
		// first number is file width and sencond one is height
		else {
			//size = str_cat(size, input);
			// width is not assigned
			if (width == 0) {
				width = atoi(read_line);
			}
			else {
				height = atoi(read_line);
			}
		}
	}

	return fp;
}

/** Read data from the file and do pre-processing
Store the pixel data into the array and return the pointer of the array
*/
unsigned char **read_data(const char *fname) {
	FILE* fp;

	fp = fopen(fname, "rb");
	if (fp == NULL) { perror(fname); return 0; }

	// read header
	fp = read_header(fp);

	if (c > width || c > height) { 
		printf("\nThe value of c is invalide"); 
		
		exit(0); 
	}

	unsigned char **pixel_data = (unsigned char **)malloc(height * sizeof(unsigned char *)); // the memory allocate to store the pixel data

	if (strcmp(ftype, "P3") == 0) {
		for (int row = 0; row < height; row++) {
			pixel_data[row] = (unsigned char *)malloc(width * 3 * sizeof(unsigned char));
		}
		unsigned char *term = (unsigned char *)malloc(sizeof(unsigned char) * 1);
		int i = 0;
		int row, col;
		while (fscanf(fp, "%u", &term) == 1) {
			row = i / (width * 3);
			col = i % (width * 3);
			(*(pixel_data + row))[col] = (unsigned char)term;
			i++;
		}
		fclose(fp);
	}

	if (strcmp(ftype, "P6") == 0) {
		int column, row, k;
		unsigned char * buf = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));

		fread(buf, sizeof(unsigned char), width*height * 3, fp); // read all data from the file
		for (row = 0, k = 0; row < height; row++) {
			pixel_data[row] = (unsigned char *)malloc(width * 3 * sizeof(unsigned char));
			for (column = 0; column < width * 3; column++, k++) {
				*(*(pixel_data + row) + column) = (unsigned int)buf[k];
			}
		}
		free(buf);
		fclose(fp);
	}


	return pixel_data;
}

inline double log2(double n) {
	return log(n) / log(2);
}

void cpu_cal() {
	printf("CPU RUNNING\n");
	int i, j, ci, cj; // for index
	int r_ = 0, g_ = 0, b_ = 0; // to calculate the average rgb
	int r_acc = 0, g_acc = 0, b_acc = 0; // accumulated rgb for each block
	int rc = 0, gc = 0, bc = 0; // accumulated rgb for whole image
	int i_c = c, j_c = c; // to solve the boundry overflow problem
	int counter;

	for (i = 0; i < height; i += c) { // row in image
		for (j = 0; j < width * 3; j += 3 * c) { // column in image

			for (ci = i, r_acc = 0, g_acc = 0, b_acc = 0, counter = 0; ci < i + c && ci < height; ci++) {  // row in block
				for (cj = j; cj < j + c * 3 && cj < width * 3; cj += 3, counter++) {  // column in block
					r_acc += *(*(data + ci) + cj + 0);
					g_acc += *(*(data + ci) + cj + 1);
					b_acc += *(*(data + ci) + cj + 2);
				}
			}
			unsigned int
				r_avg = r_acc / counter,
				g_avg = g_acc / counter,
				b_avg = b_acc / counter;

			rc += r_acc;
			gc += g_acc;
			bc += b_acc;

			for (ci = i; ci < i + c && ci < height; ci++) {  // row in block
				for (cj = j; cj < j + c * 3 && cj < width * 3; cj += 3) {  // column in block

					*(*(data + ci) + cj + 0) = r_avg;
					*(*(data + ci) + cj + 1) = g_avg;
					*(*(data + ci) + cj + 2) = b_avg;
				}
			}

			r_ += r_avg;
			g_ += g_avg;
			b_ += b_avg;

		}
	}

	r = rc / (width * height);
	g = gc / (width * height);
	b = bc / (width * height);
}

void openmp_cal() {
	printf("OPENMP RUNNING\n");

	int r_ = 0, g_ = 0, b_ = 0; // to calculate the average rgb
	int rc = 0, gc = 0, bc = 0; // accumulated rgb for whole image
	int i;
	int r_acc, g_acc, b_acc; // accumulated rgb

#pragma omp parallel for reduction(+: r_ , g_ , b_)
	for (i = 0; i < height; i += c) { // row in image
		int j;
		int  ci, cj; // for index
		int counter = 0;
		int r_avg = 0, g_avg = 0, b_avg = 0;

#pragma omp parallel for reduction(+: rc , gc , bc,r_acc, g_acc, b_acc)
		for (j = 0; j < width * 3; j += 3 * c) { // column in image

			for (ci = i, r_acc = 0, g_acc = 0, b_acc = 0, counter = 0; ci < i + c && ci < height; ci++) {  // row in block
				for (cj = j; cj < j + c * 3 && cj < width * 3; cj += 3, counter++) {  // column in block

					r_acc += *(*(data + ci) + cj + 0);
					g_acc += *(*(data + ci) + cj + 1);
					b_acc += *(*(data + ci) + cj + 2);
				}
			}
			r_avg = r_acc / counter;
			g_avg = g_acc / counter;
			b_avg = b_acc / counter;

			rc += r_acc;
			gc += g_acc;
			bc += b_acc;

			for (ci = i; ci < i + c && ci < height; ci++) {  // row in block
				for (cj = j; cj < j + c * 3 && cj < width * 3; cj += 3) {  // column in block
					*(*(data + ci) + cj + 0) = r_avg;
					*(*(data + ci) + cj + 1) = g_avg;
					*(*(data + ci) + cj + 2) = b_avg;
				}
			}
			r_ += r_avg;
			g_ += g_avg;
			b_ += b_avg;
		}
	}
	r = rc / (width * height);
	g = gc / (width * height);
	b = bc / (width * height);

}


/*Pixcel based add up value, per pixcel per thread*/
__global__
void add_up(unsigned char *data, int width, int height, int c, int new_width, int new_height, unsigned int * add_up_data, unsigned int * c_array, unsigned long long int * rgb_all) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < width*height) {
		int loc_row = i / width;
		int loc_col = i % width;
		int loc_row_new = loc_row / c;
		int loc_col_new = loc_col / c;

		atomicAdd((add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 0), *(data + i * 3 + 0));
		atomicAdd((add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 1), *(data + i * 3 + 1));
		atomicAdd((add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 2), *(data + i * 3 + 2));
		atomicAdd((c_array + (loc_row_new * new_width + loc_col_new)), 1); // to count how many pixcel in a mosic block
		atomicAdd((rgb_all + 0), *(data + i * 3 + 0)); // to addup all rgb value
		atomicAdd((rgb_all + 1), *(data + i * 3 + 1)); // to addup all rgb value 
		atomicAdd((rgb_all + 2), *(data + i * 3 + 2)); // to addup all rgb value 
	}
}

/*Mosaic based add up value, per mosic cell per block*/
__global__
void add_up_optimised(unsigned char *data, int width, int height, int c, int new_width, int new_height, unsigned int * add_up_data, unsigned int * c_array, unsigned long long int * rgb_all, int per_mosaic_block_num) {
	__shared__  unsigned int r;
	__shared__  unsigned int g;
	__shared__  unsigned int b;

	int i = (threadIdx.x / c + blockIdx.y*c)*width + (blockIdx.x*c + threadIdx.x % c);
	//					row	* width				 +					col		 	

	int MAXIMUM_WIDTH = c > 32 ? 32 : c;
	int loc_row = i / width;
	int loc_col = i % width;
	int loc_row_new = loc_row / c;
	int loc_col_new = loc_col / c;

	int blockid = blockIdx.x + gridDim.x*blockIdx.y;
	int cellid = blockIdx.x / per_mosaic_block_num + (blockIdx.y / per_mosaic_block_num)*(gridDim.x / per_mosaic_block_num);
	int capacity = c * c;

	if (per_mosaic_block_num > 1) {
		if (blockIdx.x % per_mosaic_block_num == per_mosaic_block_num - 1) {
			if (blockIdx.y % per_mosaic_block_num == per_mosaic_block_num - 1) {
				capacity = c - (per_mosaic_block_num - 1) * 32;
				capacity = capacity * capacity;
			}
			capacity = (c - (per_mosaic_block_num - 1) * 32) * 32;
		}
	}

	if (threadIdx.x < capacity-1) {
		printf("%d %d %d %d %d %d \n", i, cellid, blockid, *(data + i * 3 + 0), *(data + i * 3 + 1), *(data + i * 3 + 2));
		atomicAdd(&r, *(data + i * 3 + 0));
		atomicAdd(&g, *(data + i * 3 + 1));
		atomicAdd(&b, *(data + i * 3 + 2));	
	}	

	__syncthreads();

	if (threadIdx.x == 0) {
		printf("---%d %d %d %d %d %d \n", i, cellid, blockid, r, g, b);
		atomicAdd((add_up_data + cellid * 3 + 0), r);
		atomicAdd((add_up_data + cellid * 3 + 1), g);
		atomicAdd((add_up_data + cellid * 3 + 2), b);
		atomicAdd((c_array + cellid), capacity); // to count how many pixcel in a mosic block

		atomicAdd((rgb_all + 0), r); // to addup all rgb value
		atomicAdd((rgb_all + 1), g); // to addup all rgb value 
		atomicAdd((rgb_all + 2), b); // to addup all rgb value 
	}
}

/*calculate the average value in mosaic cell and replace the original value by the value in mosaic cell*/
__global__
void avg(unsigned char * data, int width, int height, int c, int new_width, int new_height, unsigned int * add_up_data, unsigned int * c_array) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < width*height) {
		int loc_row = i / width;
		int loc_col = i % width;
		int loc_row_new = loc_row / c;
		int loc_col_new = loc_col / c;

		*(data + i * 3 + 0) = (*(add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 0) / *(c_array + (loc_row_new * new_width + loc_col_new)));
		*(data + i * 3 + 1) = (*(add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 1) / *(c_array + (loc_row_new * new_width + loc_col_new)));
		*(data + i * 3 + 2) = (*(add_up_data + (loc_row_new * new_width + loc_col_new) * 3 + 2) / *(c_array + (loc_row_new * new_width + loc_col_new)));
	}
}

unsigned char * gpu_cal(unsigned char *gpu_data) {
	size_t size = height * width * 3 * sizeof(unsigned char);
	int add_up_data_width;
	int add_up_data_height;
	unsigned int *add_up_data_dev, *add_up_data_host; // to calculate the total rgb value in a mosic cell 
	unsigned int *c_array_dev, *c_array_host; // to count how many pixels in a mosic cell

	unsigned char *data_1d_dev; // image data
	unsigned long long int *rgb_all_dev, *rgb_all_host; // all rgb value addup
	const int BLOCK_SIZE = 512;

	add_up_data_width = width % c == 0 ? width / c : (width / c + 1);
	add_up_data_height = height % c == 0 ? height / c : (height / c + 1);

	add_up_data_host = (unsigned int *)malloc(add_up_data_width * add_up_data_height * 3 * sizeof(unsigned int));
	c_array_host = (unsigned int *)malloc(add_up_data_width * add_up_data_height * sizeof(unsigned int));
	rgb_all_host = (unsigned long long int *)malloc(3 * sizeof(unsigned long long int));

	cudaMalloc(&data_1d_dev, size);
	cudaMalloc(&add_up_data_dev, add_up_data_width * add_up_data_height * 3 * sizeof(unsigned int));
	cudaMalloc(&c_array_dev, add_up_data_width * add_up_data_height * sizeof(int));
	cudaMalloc(&rgb_all_dev, 3 * sizeof(unsigned long long int));
	cudaMemcpy(data_1d_dev, gpu_data, size, cudaMemcpyHostToDevice);


	// excutive addup kernel function
	add_up << < ((size / 3) / BLOCK_SIZE) > 0 ? (size / 3) / BLOCK_SIZE : 1, BLOCK_SIZE >> > (data_1d_dev, width, height, c, add_up_data_width, add_up_data_height, add_up_data_dev, c_array_dev, rgb_all_dev);
	
	// excutive average kernel function
	avg << < ((size / 3) / BLOCK_SIZE) > 0 ? (size / 3) / BLOCK_SIZE : 1, BLOCK_SIZE >> > (data_1d_dev, width, height, c, add_up_data_width, add_up_data_height, add_up_data_dev, c_array_dev);

	cudaMemcpy(rgb_all_host, rgb_all_dev, 3 * sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
	cudaMemcpy(gpu_data, data_1d_dev, width * height * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	r = *(rgb_all_host + 0) / width / height;
	g = *(rgb_all_host + 1) / width / height;
	b = *(rgb_all_host + 2) / width / height;

	cudaFree(add_up_data_dev);
	free(add_up_data_host);

	cudaFree(c_array_dev);
	free(c_array_host);

	cudaFree(data_1d_dev);

	cudaFree(rgb_all_dev);
	free(rgb_all_host);

	return gpu_data;
}

unsigned char * gpu_cal_optimised(unsigned char *gpu_data) {
	unsigned int *add_up_data_dev, *add_up_data_host; // to calculate the total rgb value in a mosic cell 
	unsigned int *c_array_dev, *c_array_host; // to count how many pixels in a mosic cell
	unsigned char *data_1d_dev; // image data
	unsigned long long int *rgb_all_dev, *rgb_all_host; // all rgb value addup
	const int BLOCK_SIZE = 512;

	int add_up_data_width = width % c == 0 ? width / c : (width / c + 1);
	int add_up_data_height = height % c == 0 ? height / c : (height / c + 1);
	int BLOCKBIM = c > 32 ? 32 : c; // maximun is 32
	int BLOCK_PER_MOSAIC = c / ( BLOCKBIM + 1 ) + 1;
	
	dim3 block(BLOCKBIM*BLOCKBIM, 1, 1);
	dim3 grid(BLOCK_PER_MOSAIC * add_up_data_width, BLOCK_PER_MOSAIC * add_up_data_height, 1);
	size_t size = height * width * 3 * sizeof(unsigned char);

	add_up_data_host = (unsigned int *)malloc(add_up_data_width * add_up_data_height * 3 * sizeof(unsigned int));
	c_array_host = (unsigned int *)malloc(add_up_data_width * add_up_data_height * sizeof(unsigned int));
	rgb_all_host = (unsigned long long int *)malloc(3 * sizeof(unsigned long long int));

	cudaMalloc(&data_1d_dev, size);
	cudaMalloc(&add_up_data_dev, add_up_data_width * add_up_data_height * 3 * sizeof(unsigned int));
	cudaMalloc(&c_array_dev, add_up_data_width * add_up_data_height * sizeof(int));
	cudaMalloc(&rgb_all_dev, 3 * sizeof(unsigned long long int));
	cudaMemcpy(data_1d_dev, gpu_data, size, cudaMemcpyHostToDevice);

	add_up_optimised << < grid, block >> > (data_1d_dev, width, height, c, add_up_data_width, add_up_data_height, add_up_data_dev, c_array_dev, rgb_all_dev, BLOCK_PER_MOSAIC);
	
	avg << < ((size / 3) / BLOCK_SIZE) > 0 ? (size / 3) / BLOCK_SIZE : 1, BLOCK_SIZE >> > (data_1d_dev, width, height, c, add_up_data_width, add_up_data_height, add_up_data_dev, c_array_dev);

	cudaDeviceSynchronize();

	cudaMemcpy(rgb_all_host, rgb_all_dev, 3 * sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
	cudaMemcpy(c_array_host, c_array_dev, add_up_data_width * add_up_data_height * sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(gpu_data, data_1d_dev, width * height * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	r = *(rgb_all_host + 0) / width / height;
	g = *(rgb_all_host + 1) / width / height;
	b = *(rgb_all_host + 2) / width / height;

	cudaFree(add_up_data_dev);
	free(add_up_data_host);

	cudaFree(c_array_dev);
	free(c_array_host);

	cudaFree(data_1d_dev);

	cudaFree(rgb_all_dev);
	free(rgb_all_host);

	return gpu_data;
}

int output(char * fname) {
	FILE* fp;
	int row, column, p_i, index, i;
	char* all_data;
	unsigned char* bin_data;
	char str_buf[10];
	char* char_num = (char*)malloc(4);
	int s = 0;

	printf("\nStart writing---------------\n");

	switch (output_mode) {

	case(PPM_PLAIN_TEXT):

		fp = fopen(fname, "w");
		fputs("P3\n", fp);
		fputs("# COM6521 Assignment test output\n", fp);
		sprintf(str_buf, "%d\n", width);
		fputs(str_buf, fp);
		sprintf(str_buf, "%d\n", height);
		fputs(str_buf, fp);
		sprintf(str_buf, "%d\n", 255);
		fputs(str_buf, fp);

		// format all data into string and write it into file 
		all_data = (char *)malloc(width*height * 13 * sizeof(char));
		memset(all_data, '\0', width*height * 13 * sizeof(char));

		for (row = 0, p_i = 0, index = 0; row < height; row++, p_i++, index++) {
			for (column = 0; column < width * 3; column++, i = 0, index++) {   // process number by number
				sprintf(char_num, "%d\0", *(*(data + row) + column));
				for (i = 0; *(char_num + i) != '\0' && i < 3; i++, index++) {
					*(all_data + index) = *(char_num + i);
				}
				if (p_i == 3) {
					*(all_data + index) = '\t';
					p_i = 0;
				}
				else {
					*(all_data + index) = ' ';
				}
			}
			*(all_data + index) = '\n';
		}
		fputs(all_data, fp);
		free(all_data);

		fclose(fp);
		break;


	case(PPM_BINARY):
		fp = fopen(fname, "wb");
		fputs("P6\n", fp);
		fputs("# COM6521 Assignment test output\n", fp);
		sprintf(str_buf, "%d\n", width);
		fputs(str_buf, fp);
		sprintf(str_buf, "%d\n", height);
		fputs(str_buf, fp);
		sprintf(str_buf, "%d\n", 255);
		fputs(str_buf, fp);
		bin_data = (unsigned char*)malloc(width*height * 3 * sizeof(unsigned char));

		for (row = 0, index = 0; row < height; row++) {
			for (column = 0; column < width * 3; column++, index++) {
				*(bin_data + index) = (unsigned char)*(*(data + row) + column);
			}
		}

		fwrite(bin_data, sizeof(unsigned char), width*height * 3 * sizeof(unsigned char), fp);
		fclose(fp);
		free(bin_data);

		break;
	}
	printf("The file has been saved as %s", out_file);
	return SUCCESS;
}


