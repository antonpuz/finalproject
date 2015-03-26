#include <stdio.h>
#include <jpeglib.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <setjmp.h>
#include <math.h>
#include <string.h>

/*  we will be using this uninitialized pointer later to store raw, uncompressd image */
unsigned char *raw_image = NULL;

#define SHIT 1
#define INPUT "test7.jpg"
typedef JCOEF FAR DCTCoeff;


/***** System management defines *****/

#define NOISE_MEAN gNoiseMean
int		gNoiseMean;
#define NOISE_MAX  gNoiseMax
int		gNoiseMax;

#define LUMINANCE_CHANGE gLuminance
int		gLuminance;
#define MIDDLE_LOW_FREQ_FACTOR gMiddleLowFrequencies
int		gMiddleLowFrequencies;
#define LOW_FREQ_FACTOR gLowFrequencyFactor
int		gLowFrequencyFactor;
#define LOW_PARAMS_TO_MOD gLowParamsToMod
int		gLowParamsToMod;

#define LAMBDA gLambda
double	gLambda;

#define COEF_TO_CLEAR gCorfficientsToClear
int		gCorfficientsToClear;

#define CONTRAST_DEF gContrastDef
double	gContrastDef;

#define ANOMALITY_STD_VALUE gAnomalitySTD
int		gAnomalitySTD;

/*************************************/

/* dimensions of the image we want to write */
int width = 960;
int height = 720;
int bytes_per_pixel = 3; /* or 1 for GRACYSCALE images */
int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images OR JCS_RGB for color images*/

int Seq [64] = {0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63};
int * Saved;

int coeffPlaces[][20]={
	{1,0},
	{2,1,8},
	{3,2,9,16},
	{4,3,10,17,24},
	{5,4,11,18,25,32},
	{6,5,12,19,26,33,40},
	{7,6,13,20,27,34,41,48},
	{8,7,14,21,28,35,42,49,56},
	{7,15,22,29,36,43,50,57},
	{6,23,30,37,44,51,58},
	{5,31,38,45,52,59},
	{4,39,46,53,60},
	{3,47,54,61},
	{2,55,62},
	{1,63}
};

int blocks = 0;
double stdArray[3][64];
long long meanArray[3][64];
long long meanHezka2Array[3][64];

/**
 * write_jpeg_file Writes the raw image data stored in the raw_image buffer
 * to a jpeg image with default compression and smoothing options in the file
 * specified by *filename.
 *
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to save to
 *
 */
int read_jpeg_file(char *filename) ;
void sumAll(JCOEFPTR *blockptr_one, int component);
void create_std_array();
void init(); //get the paramters from the INI file

void removeAnomalityByVariance(JCOEFPTR *blockptr_one, int component);
void removeAnomalityByVariance(JCOEFPTR *blockptr_one, int component)
{
	if(blockptr_one == NULL)return;
	int bi;
	for(bi=8 ; bi<64 ; bi++)
	{
		if(abs((*blockptr_one)[Seq[bi]] - meanArray[component][bi]) >  ANOMALITY_STD_VALUE*stdArray[component][bi] )
		{
			//printf("old value %d\n", (*blockptr_one)[Seq[bi]]);
			(*blockptr_one)[Seq[bi]] = (short)meanArray[component][bi];
			//printf("The value was: %d the std: %f changed to %d\n", (*blockptr_one)[Seq[bi]], stdArray[component][bi], (*blockptr_one)[Seq[bi]]);
		}
	}
}

void calculate_mean_var(JCOEFPTR *blockptr_one,int component);
void calculate_mean_var(JCOEFPTR *blockptr_one,int component)
{
	return;
}

void maximAnom(JCOEFPTR *blockptr_one,int component);
void maximAnom(JCOEFPTR *blockptr_one,int component)
{
	int bi;
	if(blockptr_one == NULL)return;
	if(!component)
	{
		for(bi=5 ; bi<64 ; bi++){
			if(bi>=2 && bi<=62){
				if(
				(*blockptr_one)[Seq[bi]]>((*blockptr_one)[Seq[bi-1]]+(*blockptr_one)[Seq[bi-2]])/2
				&&
				(*blockptr_one)[Seq[bi]]>((*blockptr_one)[Seq[bi+1]]+(*blockptr_one)[Seq[bi+2]])/2
				)
				// ANOMALITY CHANGE
				{
					(*blockptr_one)[Seq[bi]] -=
					(((*blockptr_one)[Seq[bi-1]]+(*blockptr_one)[Seq[bi-2]])
					+
					((*blockptr_one)[Seq[bi+1]]+(*blockptr_one)[Seq[bi+2]]))/4;

				}
			}

		}

	}

}

void clearAnomality(JCOEFPTR *blockptr_one);
void clearAnomality(JCOEFPTR *blockptr_one)
{
	int bi;
	int real_place;
	int values_accumulator;
	for(bi=11 ; bi<63 ; bi++)
	{
		real_place = Seq[bi]; // the place of the coefficient in the zig zag
		//if(((*blockptr_one)[real_place] > (*blockptr_one)[Seq[bi-1]]) && (*blockptr_one)[real_place] > (*blockptr_one)[Seq[bi+1]]) //anomality detected!!
		if((((*blockptr_one)[real_place] > (*blockptr_one)[Seq[bi-1]]) && ((*blockptr_one)[real_place] > (*blockptr_one)[Seq[bi+1]])/* && ((*blockptr_one)[real_place]>0)*/) ||
			(((*blockptr_one)[real_place] < (*blockptr_one)[Seq[bi-1]]) && ((*blockptr_one)[real_place] < (*blockptr_one)[Seq[bi+1]])/* && ((*blockptr_one)[real_place] <0)*/)) //anomality detected!!
		{
			//printf("Anomaly detected!\n");
			values_accumulator = ((*blockptr_one)[Seq[bi+1]] + (*blockptr_one)[Seq[bi-1]])/2;
			//values_accumulator = values_accumulator/2;
			//(*blockptr_one)[Seq[bi]] = values_accumulator;
			(*blockptr_one)[Seq[bi]] = (4*values_accumulator+0*(*blockptr_one)[Seq[bi]])/4;
		}
		//(*blockptr_one)[Seq[bi]] = (*blockptr_one)[Seq[bi]]*CONTRAST_DEF;
	}
}

void addContrast(JCOEFPTR *blockptr_one);
void addContrast(JCOEFPTR *blockptr_one)
{
	if(blockptr_one == NULL)return;
	int bi;
	for(bi=8 ; bi<64 ; bi++)
	{
		(*blockptr_one)[Seq[bi]] = (*blockptr_one)[Seq[bi]]*CONTRAST_DEF;
	}
}


void clearDCTCoeffs(JCOEFPTR *blockptr_one);
void clearDCTCoeffs(JCOEFPTR *blockptr_one)
{
int n,t;
for (n=0 ; n < COEF_TO_CLEAR ; n++)
    Saved[n]=-1;
for (n=0 ; n < COEF_TO_CLEAR ; n++)
{
    do
    {
        t= Seq[10+rand()%54];
    }while(NewItem(t)==1);
    Saved[n]=t;
}
if(blockptr_one == NULL)return;
    int bi;


    for(bi=0 ; bi<COEF_TO_CLEAR ; bi++)
    {
        (*blockptr_one)[Saved[bi]] = 0;
    }

}

void quantization(JCOEFPTR *blockptr_one);
void quantization(JCOEFPTR *blockptr_one)
{
if(blockptr_one == NULL)return;
int STARTING_E = 4;
int V=0;
int P=0;
int E=0;
//This is the quantization table that will be used to quantize the image
int QTable[15] = {1,1,1,1,3,3,6,6,6,10,10,15,15,15,15};
// THERE ARE 15 E's
    for(E=STARTING_E;E<15;E++){
        V = coeffPlaces[E][0];
        for(P=1;P<(V+1);P++){
            (*blockptr_one)[coeffPlaces[E][P]] = (*blockptr_one)[coeffPlaces[E][P]]/QTable[E];
            (*blockptr_one)[coeffPlaces[E][P]] = (*blockptr_one)[coeffPlaces[E][P]]*QTable[E];
        }

    }

}

void enchanceImage(JCOEFPTR *blockptr_one);
void enchanceImage(JCOEFPTR *blockptr_one)
{
	//calculate the original En
	//DCTCoeff dctArray[15] = {0};
	double enchanceddctArray[15] = {0};
	double dctArray[15] = {0};
	double Hfactor[15] = {0};
	JCOEF tmpDataArray[64] = {0};
	//for loop
	int index, innerIndex;
	int dctShift = 0;
	for(index=0 ; index<64 ; index+=8)
	{
		for(innerIndex = 0 ; innerIndex < 8 ; innerIndex++)
		{
			dctArray[dctShift+innerIndex] += abs((*blockptr_one)[index+innerIndex]);
			//printf("%d ", (*blockptr_one)[index+innerIndex]);
		}
		dctShift++;
	}
	//printf("\n");
	for(index = 0 ; index < 15 ; index++)
	{
		//devide the summed coefficient values by the number of numbers added
		if(index < 8)
		{
			dctArray[index] /= (double)(index+1);
		}
		else
		{
			dctArray[index] /= (double)(14-index+1);
		}
	}
	//calculate the enhanced En	
	enchanceddctArray[0] = dctArray[0];
	Hfactor[0] = 1;
	double regularSum;
	double enchancedSum;
	for(index=1 ; index<15 ; index++)
	{
		regularSum = 0;
		enchancedSum = 0;
		for(innerIndex = 0 ; innerIndex < index ; innerIndex++)
		{

			regularSum += dctArray[innerIndex];
			enchancedSum += enchanceddctArray[innerIndex];
			//printf("adding reg %f so suym is: %f, enchanced added: %f and sum: %f\n",dctArray[innerIndex], regularSum, enchanceddctArray[innerIndex], enchancedSum);
		}
		if(regularSum != 0)
			Hfactor[index] = enchancedSum/regularSum;
		else
			Hfactor[index] = 1;
		//printf("reg array val: %f encahanced: %f and enchanced: %f\n", regularSum, enchancedSum, Hfactor[index]);
		enchanceddctArray[index] = LAMBDA*Hfactor[index]*dctArray[index];
	}
	//change the coeffs
	dctShift = 0;
//	printf("new coeff\n");
	for(index=0 ; index<64 ; index+=8)
	{
		for(innerIndex = 0 ; innerIndex < 8 ; innerIndex++)
		{
			if(index == 0 && innerIndex == 0)
			{
				//printf("new coeff %d \n", (*blockptr_one)[index+innerIndex]);
				tmpDataArray[0] = (*blockptr_one)[0];
				continue;
			}
			//printf("calculating index: %d, inner %d, dctshift %d\n",index,innerIndex,dctShift);
			//printf("old value: %d lambda:%f enhance:%f \n", (*blockptr_one)[index+innerIndex], lambda, Hfactor[dctShift + innerIndex]);
			int will_wrt;
			double tmp = LAMBDA*Hfactor[dctShift + innerIndex]*((*blockptr_one)[index+innerIndex]);
			if(tmp > 32750) tmp = 32750;
			if(tmp < -32750) tmp = -32750;
			//will_wrt = (int)tmp;
			//printf("writing data in cell %d\n", index+innerIndex);
			(*blockptr_one)[index+innerIndex] = (short)tmp;
			//(*blockptr_one)[index+innerIndex] = (short)tmp;
			//printf("%d(%f) ", tmpDataArray[index+innerIndex], tmp);
		}
		dctShift++;
	}
//	for(index=0 ; index<64 ; index++)
//	{
//		//if(index < 10)
//		(*blockptr_one)[index] = (short)tmpDataArray[index];
//		//(*blockptr_one)[index] = (*blockptr_one)[index];
//		printf("%d ", tmpDataArray[index]);
//		if(tmpDataArray[index] != (*blockptr_one)[index])
//		{
//			printf("############%d############ ", tmpDataArray[index]);
//			exit(1);
//		}
//
//	}
//	printf("\n");

}

void rotateDCTCoeff(JCOEFPTR *blockptr_one);
void rotateDCTCoeff(JCOEFPTR *blockptr_one)
{
	if(blockptr_one == NULL)return;
	int index;
	int tmp;
	int coeffTmp;
	int rnd;
	int j;
	//coeffPlaces[15][];
//	printf("The old values: ");
//	for(j = 0 ; j < 64 ; j++)
//	{
//		printf("%d, ", (*blockptr_one)[j]);
//	}
//	printf("\n");
	for(index = 4 ; index < 15 ; index++)
	{
		//rnd = rand() % coeffPlaces[index][0]; // rotate according to the number of coeffs in the block

		j = coeffPlaces[index][0];
		while(j > 0)
		{
			rnd = rand() % j;
			tmp = coeffPlaces[index][1 + rnd];

			coeffTmp = (*blockptr_one)[coeffPlaces[index][1 + rnd]];
			//printf("1: %d, 2: %d\n", tmp, coeffTmp);
			(*blockptr_one)[coeffPlaces[index][1 + rnd]] = (*blockptr_one)[coeffPlaces[index][1]];
			(*blockptr_one)[coeffPlaces[index][1]] = coeffTmp;

			coeffPlaces[index][1 + rnd] = coeffPlaces[index][1];
			coeffPlaces[index][1] = tmp;
			j--;
		}

	}
//	printf("The new values: ");
//	for(j = 0 ; j < 64 ; j++)
//	{
//		printf("%d, ", (*blockptr_one)[j]);
//	}
//	printf("\n");
}

void addNoiseDCT(JCOEFPTR *blockptr_one, int component);
void addNoiseDCT(JCOEFPTR *blockptr_one, int component)
{
	if(blockptr_one == NULL)return;
	int noise_var;
	int bi;
	int ra;
	int noise;
	for(bi=0 ; bi<64 ; bi++)
	{
		if(bi >= 15){
			//printf("cal noise var\n");
			noise_var = (!component) ? NOISE_MEAN + (rand() % (2*NOISE_MAX)) - NOISE_MAX : NOISE_MEAN + (rand() % (2*((int)ceil((double)NOISE_MAX/(double)3)))) - (NOISE_MAX/3);
			//printf("finish noise var calc\n");
			//noise = NOISE_MEAN + (rand() % (2*NOISE_MAX)) - NOISE_MAX;
			(*blockptr_one)[Seq[bi]] = (*blockptr_one)[Seq[bi]]+noise_var;
		}
	}
}

void printDCTCoeff(JCOEFPTR *blockptr_one);
void printDCTCoeff(JCOEFPTR *blockptr_one)
{
	if(blockptr_one == NULL)return;

	int bi;
	printf("Print coeffs\n");
	for(bi=0 ; bi<64 ; bi++)
	{
		printf("%d ", (*blockptr_one)[bi]);
	}
	printf("\n");
}

void changeLuminance(JCOEFPTR *blockptr_one, int component);
void changeLuminance(JCOEFPTR *blockptr_one, int component)
{
	if(blockptr_one == NULL)return;

	int i, rnd;

	int very_low[] = {1,2,8,9,16};
	int middle_low[] = {3,10,17,24};

	if(!component)
	{
//		printf("The old DCT val: %d\n", (*blockptr_one)[0]);
		(*blockptr_one)[0] = (*blockptr_one)[0] + LUMINANCE_CHANGE;
//		printf("The new DCT val: %d\n", (*blockptr_one)[0]);
	}


	for (i = 0 ; i < (sizeof(middle_low)/sizeof(int)) ; i++)
	{
		(*blockptr_one)[middle_low[i]] = (*blockptr_one)[middle_low[i]] *= MIDDLE_LOW_FREQ_FACTOR;
	}

	for (i = 0 ; i < LOW_PARAMS_TO_MOD ; i++)
	{
		rnd = rand() % (sizeof(very_low)/sizeof(int));
		(*blockptr_one)[rnd] = (*blockptr_one)[rnd] *= LOW_FREQ_FACTOR;
	}

}


int NewItem();
int NewItem(int t){
int k;
    for(k=0;k<COEF_TO_CLEAR;k++)
        if(Saved[k]==t)
            return 1;
return 0;
}

int main(int argc, char **argv ) {
	if(argc != 3)
	{
		printf("wrong number of params\n");
		exit(0);
	}

	init();
	//exit(1);

	memset(meanArray, 0, 3*64*sizeof(long long));
	memset(meanHezka2Array, 0, 3*64*sizeof(long long));
	memset(stdArray, 0, 3*64*sizeof(double));

	//char *infilename = INPUT, *outfilename = "test_out.jpg";
	char *infilename = argv[1], *outfilename = argv[2];
	//read_jpeg_file(infilename);
	/* Try opening a jpeg*/
//	if (read_jpeg_file(infilename) > 0) {
//		/* then copy it to another file */
//		if (write_jpeg_file(outfilename) < 0)
//			return -1;
//	} else
//		return -1;
		/* yup, we succeeded! */
	srand (time(NULL));

	read_coefficients(infilename,outfilename);

	return 0;
}


int read_coefficients(char* filename,char* output)
{

	  /* This struct contains the JPEG decompression parameters and pointers to
	   * working space (which is allocated as needed by the JPEG library).
	   */
	  struct jpeg_decompress_struct cinfo;

	  struct jpeg_error_mgr jerr;
	  /* More stuff */
	  FILE * infile;        /* source file */

	  /* In this example we want to open the input file before doing anything else,
	   * so that the setjmp() error recovery below can assume the file is open.
	   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
	   * requires it in order to read binary files.
	   */

	  if ((infile = fopen(filename, "rb")) == NULL) {
	    fprintf(stderr, "can't open %s\n", filename);
	    return 0;
	  }

	  /* Step 1: allocate and initialize JPEG decompression object */
	  cinfo.err = jpeg_std_error(&jerr);
	  /* Now we can initialize the JPEG decompression object. */
	  jpeg_create_decompress(&cinfo);

	  /* Step 2: specify data source (eg, a file) */
	  jpeg_stdio_src(&cinfo, infile);

	  /* Step 3: read file parameters with jpeg_read_header() */
	  (void) jpeg_read_header(&cinfo, TRUE);

	  /* We can ignore the return value from jpeg_read_header since
	   *   (a) suspension is not possible with the stdio data source, and
	   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	   * See libjpeg.txt for more info.
	   */

	//Createa a compress file that will be used to write changed DCT to it - same steps as for the decompress
	struct jpeg_compress_struct cinfo2;
	struct jpeg_error_mgr jerr2;
	jpeg_create_compress(&cinfo2);

	/* this is a pointer to one row of image data */
	JSAMPROW row_pointer2[1];
	FILE *outfile = fopen(output, "wb");

	if (!outfile) {
		printf("Error opening output jpeg file %s\n!", "bananas.jpeg");
		return -1;
	}
	cinfo2.err = jpeg_std_error(&jerr2);
	jpeg_create_compress(&cinfo2);
	jpeg_stdio_dest(&cinfo2, outfile);

	//set the number of components for the default params
	cinfo2.input_components = bytes_per_pixel;

	/* default compression parameters, we shouldn't be worried about these */
	jpeg_set_defaults(&cinfo2); //TODO: CHeck if this is necessary with the read coefficients function called later.

	/* Setting the parameters of the output file here */
	/* We use the quality of the first image on the second*/
	/* jpeg_set_quality(&cinfo2, 50, TRUE); */
	//cinfo2.dct_method = JDCT_FLOAT;
	//cinfo2.num_components = 3;
	//cinfo2.in_color_space = color_space;
	//cinfo2.image_width = width;
	//cinfo2.image_height = height;
	//cinfo.data_precision = 4;

	jvirt_barray_ptr *coeffs = jpeg_read_coefficients(&cinfo);
	jpeg_copy_critical_parameters(&cinfo, &cinfo2);

	int i;
	int Changed=0;
	int TotalCoef=0;
	short TEMP[64];

	//decompress and call DCT manipulation function
	int done = 0;
	int ci;
	int by;
	int bx;
	int bi;

	for (ci = 0; ci < 3/*3*/; ci++)
	{
	    JBLOCKARRAY buffer_one;
	    JCOEFPTR blockptr_one;
	    jpeg_component_info* compptr_one;

	    compptr_one = cinfo.comp_info + ci;
	    for (by = 0; by < compptr_one->height_in_blocks; by++)
	    {
	        buffer_one = (cinfo.mem->access_virt_barray)((j_common_ptr)&cinfo, coeffs[ci], by, (JDIMENSION)1, FALSE);
	        for (bx = 0; bx < compptr_one->width_in_blocks; bx++)
	        {
	            blockptr_one = buffer_one[0][bx];

	            sumAll(&blockptr_one, ci);
	            //Maxims_function(&blockptr_one);



	        }
	    }
	}

	create_std_array();

	for (ci = 0; ci < 3/*3*/; ci++)
	{
	    JBLOCKARRAY buffer_one;
	    JCOEFPTR blockptr_one;
	    jpeg_component_info* compptr_one;

	    compptr_one = cinfo.comp_info + ci;
	    printf("height:%u, width:%u\n", compptr_one->height_in_blocks,compptr_one->width_in_blocks);

	    for (by = 0; by < compptr_one->height_in_blocks; by++)
	    {
	        buffer_one = (cinfo.mem->access_virt_barray)((j_common_ptr)&cinfo, coeffs[ci], by, (JDIMENSION)1, FALSE);

	        for (bx = 0; bx < compptr_one->width_in_blocks; bx++)
	        {
	            blockptr_one = buffer_one[0][bx];
//
	            for(i = 0 ; i < 64 ; i++)
	                TEMP[i] = blockptr_one[i];

	        	//printDCTCoeff(&blockptr_one);


	            removeAnomalityByVariance(&blockptr_one, ci);

//	        	//clearAnomality(&blockptr_one);
	            enchanceImage(&blockptr_one);
	            addContrast(&blockptr_one);
	            changeLuminance(&blockptr_one, ci);
//
//
//
//
//
	            addNoiseDCT(&blockptr_one, ci);
	        	quantization(&blockptr_one);
	        	clearDCTCoeffs(&blockptr_one);
//	        	//rotateDCTCoeff(&blockptr_one);
//
//	        	//printDCTCoeff(&blockptr_one);
//
//	        	//manipulateDCTCoeff(&blockptr_one);
	        	for(i = 0 ; i < 64 ; i++){
	        		TotalCoef++;
	        		if(TEMP[i] != blockptr_one[i])
	        			Changed++;
	        	}
	        }
	    }
	}

	printf("TotalCoef: %d\nChanged: %d\n",TotalCoef,Changed);
	double pampam;
	pampam = (double)Changed/(double)TotalCoef;
	pampam *= 100;
	printf("Percentage: %f\n",pampam);

	printf("writing the new image\n");


	jpeg_write_coefficients(&cinfo2, coeffs);
	jpeg_finish_compress(&cinfo2);
	jpeg_destroy_compress(&cinfo2);

	fclose(outfile);

	jpeg_destroy_decompress(&cinfo);

	return 1;
}

/**
 * read_jpeg_file Reads from a jpeg file on disk specified by filename and saves into the
 * raw_image buffer in an uncompressed format.
 *
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to read from
 *
 */

int read_jpeg_file(char *filename) {
	/* these are standard libjpeg structures for reading(decompression) */
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	char tmpArr[3];
	double midTmp;
	char zro = 0;
	char mid;
	int TotalSize = 151353 + 12 + 14;
	int tmp;
	short shrtTmp;
	/* libjpeg data structure for storing one row, that is, scanline of an image */
	JSAMPROW row_pointer[1];
	FILE *infile = fopen(filename, "rb");

	unsigned long location = 0;
	int i = 0;

	if (!infile) {
		printf("Error opening jpeg file %s\n!", filename);
		return -1;
	}
	/* here we set up the standard libjpeg error handler */
	cinfo.err = jpeg_std_error(&jerr);
	/* setup decompression process and source, then read JPEG header */
	jpeg_create_decompress(&cinfo);
	/* this makes the library read from infile */
	jpeg_stdio_src(&cinfo, infile);
	/* reading the image header which contains image information */
	jpeg_read_header(&cinfo, TRUE);
	/* Uncomment the following to output image information, if needed. */

	width = cinfo.image_width;
	height = cinfo.image_height;

	 printf( "JPEG File Information: \n" );
	 printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
	 printf( "Color components per pixel: %d.\n", cinfo.num_components );
	 printf( "Color space: %d.\n", cinfo.jpeg_color_space );
	 printf( "densoty unit: %d x=%u y=%d transform=%d.\n", cinfo.density_unit , cinfo.X_density, cinfo.Y_density, cinfo.Adobe_transform);

	/* Start decompression jpeg here */
	jpeg_start_decompress(&cinfo);

	/* allocate memory to hold the uncompressed image */
	raw_image = (unsigned char*) malloc(
			cinfo.output_width * cinfo.output_height * cinfo.num_components);
	/* now actually read the jpeg into the raw buffer */
	row_pointer[0] = (unsigned char *) malloc(
			cinfo.output_width * cinfo.num_components);
	/* read one scan line at a time */
	while (cinfo.output_scanline < cinfo.image_height) {
		jpeg_read_scanlines(&cinfo, row_pointer, 1);
		for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
		{
			raw_image[location++] = row_pointer[0][i];
		}

	}
	printf("The location size: %lu\n", location);
	/* wrap up decompression, destroy objects, free pointers and close open files */
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	free(row_pointer[0]);
	fclose(infile);
	/* yup, we succeeded! */
	return 1;
}

int shit()
{
	char *infilename = "test.jpg", *outfilename = "test_out.jpg";
	struct jpeg_decompress_struct cinfo;
		struct jpeg_error_mgr jerr;
		/* libjpeg data structure for storing one row, that is, scanline of an image */
		JSAMPROW row_pointer[1];

		FILE *infile = fopen(infilename, "rb");
		unsigned long location = 0;
		int i = 0;

		if (!infile) {
			printf("Error opening jpeg file %s\n!", infilename);
			return -1;
		}

		/* here we set up the standard libjpeg error handler */
		cinfo.err = jpeg_std_error(&jerr);
		/* setup decompression process and source, then read JPEG header */
		jpeg_create_decompress(&cinfo);
		/* this makes the library read from infile */
		jpeg_stdio_src(&cinfo, infile);
		/* reading the image header which contains image information */

		jpeg_read_header(&cinfo, TRUE);
		/* Uncomment the following to output image information, if needed. */

		 printf( "JPEG File Information: \n" );
		 printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
		 printf( "Color components per pixel: %d.\n", cinfo.num_components );
		 printf( "Color space: %d.\n", cinfo.jpeg_color_space );

		/* wrap up decompression, destroy objects, free pointers and close open files */
		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		free(row_pointer[0]);
		fclose(infile);

		return SHIT;
}

int write_jpeg_file(char *filename) {
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;

	/* this is a pointer to one row of image data */
	JSAMPROW row_pointer[1];
	FILE *outfile = fopen(filename, "wb");

	if (!outfile) {
		printf("Error opening output jpeg file %s\n!", filename);
		return -1;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

	/* Setting the parameters of the output file here */
	cinfo.image_width = width;
	cinfo.image_height = height;
	cinfo.input_components = bytes_per_pixel;
	cinfo.in_color_space = color_space;
	/* default compression parameters, we shouldn't be worried about these */

	jpeg_set_defaults(&cinfo);
	cinfo.num_components = 3;
	//cinfo.data_precision = 4;
	cinfo.dct_method = JDCT_FLOAT;
	jpeg_set_quality(&cinfo, 75, TRUE);
	/* Now do the compression .. */
	jpeg_start_compress(&cinfo, TRUE);
	/* like reading a file, this time write one row at a time */
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &raw_image[cinfo.next_scanline * cinfo.image_width
				* cinfo.input_components];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	/* similar to read file, clean up after we're done compressing */
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	fclose(outfile);
	/* success code is 1! */
	return 1;
}

void create_std_array()
{
	int i;

	printf("The mean array:\n");
	for(i=0;i<64;i++)
	{
//		meanHezka2Array[0][i] = round((double)meanHezka2Array[0][i]/(double)blocks);
//		meanHezka2Array[1][i] = round((double)meanHezka2Array[1][i]/(double)blocks);
//		meanHezka2Array[2][i] = round((double)meanHezka2Array[2][i]/(double)blocks);
//
//		meanArray[0][i] = round((double)meanArray[0][i]/(double)blocks);
//		meanArray[1][i] = round((double)meanArray[1][i]/(double)blocks);
//		meanArray[2][i] = round((double)meanArray[2][i]/(double)blocks);


		printf(" %lld", meanArray[0][i]);
	}
	printf("\n");

	printf("The STD array:\n");

	for(i=0;i<64;i++)
	{
//		printf("The old mean %lld\n", meanArray[0][i]);
//		printf("The mean divided %f\n", (double)meanArray[0][i]/(double)blocks);
//		printf("The old mean^2 %lld\n", meanHezka2Array[0][i]);
//		printf("The mean divided^2 %f\n", (double)meanHezka2Array[0][i]/(double)blocks);

		stdArray[0][i] = sqrt(((double)meanHezka2Array[0][i]/(double)blocks) - ((((double)meanArray[0][i])*((double)meanArray[0][i]/(double)blocks))));
		stdArray[1][i] = sqrt(((double)meanHezka2Array[1][i]/(double)blocks) - (((double)meanArray[1][i]/(double)blocks)*((double)meanArray[1][i]/(double)blocks)));
		stdArray[2][i] = sqrt(((double)meanHezka2Array[2][i]/(double)blocks) - (((double)meanArray[2][i]/(double)blocks)*((double)meanArray[2][i]/(double)blocks)));

		meanArray[0][i] = round((double)meanArray[0][i]/(double)blocks);
		meanArray[1][i] = round((double)meanArray[1][i]/(double)blocks);
		meanArray[2][i] = round((double)meanArray[2][i]/(double)blocks);
		//printf("The new mean %lld blocks %d\n", meanArray[0][i], blocks);

//
//		printf("The std %f\n", stdArray[0][i]);
	}
	//exit(1);
}


void sumAll(JCOEFPTR *blockptr_one, int component)
{

	if(blockptr_one == NULL)return;
	int bi;
	for(bi=0 ; bi<64 ; bi++)
	{
		meanArray[component][bi] += (*blockptr_one)[Seq[bi]];
		meanHezka2Array[component][bi] += ((*blockptr_one)[Seq[bi]] * (*blockptr_one)[Seq[bi]]);
	}
	if(!component)blocks++;
}

void init() //get the paramters from the INI
{
	FILE * pINI;
	char * lineBuffer = NULL;
	char * needle	= NULL;
	size_t bufferLength;
	pINI = fopen("parameters.ini", "r");
	if (pINI == NULL) {
		fprintf(stderr, "can't open parameters.ini\n");
	    return;
	}

//	while(getline(&lineBuffer, &bufferLength, pINI) > 0)
//	{
//		needle = strrchr(lineBuffer, '\n');
//		if(needle)*needle=0; //delete the new line from the buffer
//		printf("Got %s from the file with length %d\n", lineBuffer, bufferLength);
//	}
//
	/*
		#define NOISE_MEAN gNoiseMean
		int		gNoiseMean;
		#define NOISE_MAX  5

		#define LUMINANCE_CHANGE 20
		#define MIDDLE_LOW_FREQ_FACTOR 1
		#define LOW_FREQ_FACTOR 1
		#define LOW_PARAMS_TO_MOD 1

		#define LAMBDA gLambda
		double	gLambda;

		#define COEF_TO_CLEAR 30

		#define CONTRAST_DEF 1.3

		#define ANOMALITY_STD_VALUE 3

	 */

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get Noise mean
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		NOISE_MEAN = atoi(needle + 1);
		printf("The noise mean is: %d\n", NOISE_MEAN);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get Noise max
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		NOISE_MAX = atoi(needle + 1);
		printf("The noise mean is: %d\n", NOISE_MAX);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get Luminance
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		LUMINANCE_CHANGE = atoi(needle + 1);
		printf("The noise mean is: %d\n", LUMINANCE_CHANGE);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get MIDDLE_LOW_FREQ_FACTOR
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		MIDDLE_LOW_FREQ_FACTOR = atoi(needle + 1);
		printf("The noise mean is: %d\n", MIDDLE_LOW_FREQ_FACTOR);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get LOW_FREQ_FACTOR
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		LOW_FREQ_FACTOR = atoi(needle + 1);
		printf("The noise mean is: %d\n", LOW_FREQ_FACTOR);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get LOW_PARAMS_TO_MOD
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		LOW_PARAMS_TO_MOD = atoi(needle + 1);
		printf("The noise mean is: %d\n", LOW_PARAMS_TO_MOD);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get COEF_TO_CLEAR
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		COEF_TO_CLEAR = atoi(needle + 1);
		Saved = malloc(COEF_TO_CLEAR * sizeof(int));
		printf("The noise mean is: %d\n", COEF_TO_CLEAR);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get CONTRAST_DEF
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		CONTRAST_DEF = atoi(needle + 1);
		needle = strrchr(lineBuffer, '.');
		CONTRAST_DEF += (double)atoi(needle + 1)/(double)10;
		printf("The noise mean is: %f\n", CONTRAST_DEF);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get ANOMALITY_STD_VALUE
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		ANOMALITY_STD_VALUE = atoi(needle + 1);
		printf("The noise mean is: %d\n", ANOMALITY_STD_VALUE);
	}

	if(getline(&lineBuffer, &bufferLength, pINI) > 0) //get Lambda
	{
		needle = strrchr(lineBuffer, '\n');
		if(needle)*needle=0; //delete the new line from the buffer

		needle = strrchr(lineBuffer, ':');
		LAMBDA = atoi(needle + 1);
		needle = strrchr(lineBuffer, '.');
		LAMBDA += (double)atoi(needle + 1)/(double)10;
		printf("The noise mean is: %f\n", LAMBDA);
	}

	printf("The lengfth %d\n", bufferLength);


	fclose(pINI);
}
