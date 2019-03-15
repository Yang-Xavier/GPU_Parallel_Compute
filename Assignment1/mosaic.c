#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define FAILURE 0
#define SUCCESS !FAILURE

#define USER_NAME "acr18by"		//replace with your user name



typedef enum MODE { CPU, OPENMP, CUDA, ALL } MODE;
typedef enum OUTPUT_MODE { PPM_BINARY, PPM_PLAIN_TEXT } OUTPUT_MODE;

void print_help();
int process_command_line(int argc, char *argv[]);
int read_data(char* fname);
void cpu_cal();
void openmp_cal();
double log2(double n);
int output(char * fname);

unsigned int c = 0;
unsigned int width = 0;
unsigned int height = 0;
unsigned int ** data;
char* in_file;
char* out_file;
int r, g, b;
MODE execution_mode = CPU;
OUTPUT_MODE output_mode = PPM_PLAIN_TEXT;

int main(int argc, char *argv[]) {

	if (process_command_line(argc, argv) == FAILURE)
		return 1;


	//TODO: read input image file (either binary or plain text PPM) 
	data = read_data(in_file);

	//TODO: execute the mosaic filter based on the mode
	switch (execution_mode) {
	case (CPU): {
		//TODO: starting timing here
		clock_t start = clock(), diff;
		//TODO: calculate the average colour value

		cpu_cal();

		// Output the average colour value for the image
		printf("CPU Average image colour red = %d, green = %d, blue = %d \n", r, g, b);
		//TODO: end timing here
		diff = clock() - start;
		int msec = diff * 1000 / CLOCKS_PER_SEC;
		printf("CPU mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);


		break;
	}
	case (OPENMP): {
		//TODO: starting timing here
		clock_t start = clock(), diff;
		//TODO: calculate the average colour value
		openmp_cal();
		// Output the average colour value for the image
		printf("OPENMP Average image colour red = %d, green = %d, blue = %d \n", r, g, b);

		//TODO: end timing here
		diff = clock() - start;
		int msec = diff * 1000 / CLOCKS_PER_SEC;
		printf("OPENMP mode execution time took %d s and %dms\n", msec / 1000, msec % 1000);
		break;
	}
	case (CUDA): {
		printf("CUDA Implementation not required for assignment part 1\n");
		break;
	}
	case (ALL): {
		//TODO
		break;
	}


	}
	output(out_file); // Output to file

	//save the output image file (from last executed mode)
	free(data);
	getchar();
	return 0;
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

int process_command_line(int argc, char *argv[]) {
	if (argc < 7) {
		fprintf(stderr, "Error: Missing program arguments. Correct usage is...\n");
		print_help();
		return FAILURE;
	}
	//first argument is always the executable name

	//read in the non optional command line arguments
	c = (unsigned int)atoi(argv[1]);

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
	if (!strcmp(argv[8], "PPM_BINARY")) { output_mode = PPM_BINARY; };
	if (!strcmp(argv[8], "PPM_PLAIN_TEXT")) { output_mode = PPM_PLAIN_TEXT; };

	return SUCCESS;
}

/** Read data from the file and do pre-processing
Store the pixel data into the array and return the pointer of the array
*/
int read_data(const char* fname) {
	FILE* fp;
	char ppm_header[256];
	int detail[3];
	char ftype[4];
	char* rp;
	int i = 0, j = 0;
	int header_len = 0;

	fp = fopen(fname, "r");
	if (fp == NULL) { perror(fname); return 0; }

	// Read first four line and alloc the memory to the next part of data
	fgets(ftype, sizeof ftype, fp); // Get the type of the input file
	rp = strtok(ftype, "\n");	//split by using \n
	*ftype = *rp;
	header_len += strlen(ftype) + 1;

	while (j < 3) {
		char *s = fgets(ppm_header, sizeof ppm_header, fp);
		if (s[0] >= '0' && s[0] <= '9') {
			detail[j] = atoi(s);
			j++;
		}
		header_len += strlen(s);
	}

	width = (unsigned int)detail[0];
	height = (unsigned int)detail[1];
	if (c > width) { c = 1; }
	if (c > height) { c = 1; }
	unsigned int *(*pixel_data) = (unsigned int *)malloc(height * sizeof(unsigned int *)); // the memory allocate to store the pixel data


	if (strcmp(ftype, "P3") == 0) {

		int ch;
		int i = 0;
		int app_row, app_col;
		char *all = (char*)malloc(width * 13 * height);

		while ((ch = fgetc(fp)) != EOF) {  // read the data from the file and replace \t and \n with space
			if ((char)ch == '\t') {
				ch = (int)'\ ';
			}
			if ((char)ch == '\n') {
				ch = (int)'\ ';
			}
			all[i] = (char)ch;
			i++;
		}
		char* term;

		term = strtok(all, "\ "); // split the data by space and store it
		i = 0;
		int row = -1, column = 0;
		while (i < width * height * 3) {
			if (i % (width * 3) == 0) {
				row += 1;
				column = 0;
				pixel_data[row] = (unsigned int *)malloc(width * 3 * sizeof(unsigned int));
			}
			*(*(pixel_data + row) + column) = (unsigned int)atoi(term);
			term = strtok(NULL, "\ ");
			column++;
			i++;
		}
	}
	fclose(fp);

	if (strcmp(ftype, "P6") == 0) {
		fopen(fname, "rb");
		fseek(fp, header_len * sizeof(char), 0);
		// binary format
		unsigned char * buf = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
		int column, row, k;
		fread(buf, sizeof(unsigned char), width*height * 3, fp); // read all data from the file
		for (row = 0, k = 0; row < height; row++) {
			pixel_data[row] = (unsigned int *)malloc(width * 3 * sizeof(unsigned int));
			for (column = 0; column < width * 3; column++, k++) {
				*(*(pixel_data + row) + column) = (unsigned int)buf[k];
			}
		}
		free(buf);
		fclose(fp);
	}


	return pixel_data;
}

double log2(double n) {
	return log(n) / log(2);
}


void cpu_cal() {
	printf("CPU RUNNING\n");
	int i, j, ci, cj; // for index
	int r_ = 0, g_ = 0, b_ = 0; // to calculate the average rgb
	int r_acc = 0, g_acc = 0, b_acc = 0; // accumulated rgb
	int i_c = c, j_c = c; // to solve the boundry overflow problem
	int counter;
	for (i = 0; i < height; i += c) { // row in image
		for (j = 0; j < width * 3; j += 3 * c) { // column in image
												 //i_c = (i + c)  > height ? (height - c) : c; // to judge whether overflow row num
												 //j_c = (i + c * 3) > width * 3 ? (width - c) : c; // judge whether overflow column

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

	r = r_ / ((width / c)*(width / c));
	g = g_ / ((width / c)*(width / c));
	b = b_ / ((width / c)*(width / c));
}


void openmp_cal() {
	printf("OPENMP RUNNING\n");

	int r_ = 0, g_ = 0, b_ = 0; // to calculate the average rgb

	int i;
#pragma omp parallel for
	for (i = 0; i < height; i += c) { // row in image
		
		int j;
#pragma omp parallel for
		for (j = 0; j < width * 3; j += 3 * c) { // column in image
												 //i_c = (i + c)  > height ? (height - c) : c; // to judge whether overflow row num
												 //j_c = (i + c * 3) > width * 3 ? (width - c) : c; // judge whether overflow column
			int  ci, cj; // for index
			int r_acc, g_acc, b_acc; // accumulated rgb
			int i_c = c, j_c = c; // to solve the boundry overflow problem
			int counter;
			unsigned int r_avg = 0, g_avg = 0, b_avg = 0;

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
	r = r_ / ((width / c)*(width / c));
	g = g_ / ((width / c)*(width / c));
	b = b_ / ((width / c)*(width / c));

}

int output(char * fname) {
	FILE* fp;
	int row, column, p_i, index, i;
	char* all_data;
	unsigned char* bin_data;
	char* str_buf[10];
	char* char_num = (char*)malloc(4);
	int l;
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
		int s = 0;
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
		break;
	}
	printf("The file has been saved as %s", out_file);
}







