#include "DCT.hpp"

//  ---------------------------------------------------------------
// |
// |  Constants and Helper
// |
//  ---------------------------------------------------------------

#define N 8

inline float getConstantC(unsigned char x, unsigned char y) { // (2.0F / N) * getC(i) * getC(j)
	if (x != 0 && y != 0)
		// most common clause (49 / 64)
		return 0.25F;
	else if (x != 0 || y != 0)
		// second most common clause (14 / 64)
		return 0.176776695296636893184327732342353556305170059204101562F;
	else
		// least common clause (1 / 64)
		return 0.125F;
}

// float %1.25f  double %1.54f
static const float cos_2x1iPi_2N[8][8] = { // cos( (2x + 1)iπ / 2N )
	{ 1,  0.980785280403230449126182236134239036973933730893336095003F,  0.923879532511286756128183189396788286822416625863642486116F,  0.831469612302545237078788377617905756738560811987249963448F,  0.707106781186547524400844362104849039284835937688474036591F,  0.555570233019602224742830813948532874374937190754804045928F,  0.38268343236508977172845998403039886676134456248562704144F,   0.19509032201612826784828486847702224092769161775195480776F },
	{ 1,  0.831469612302545237078788377617905756738560811987249963447F,  0.38268343236508977172845998403039886676134456248562704144F,  -0.19509032201612826784828486847702224092769161775195480775F,  -0.70710678118654752440084436210484903928483593768847403659F,  -0.980785280403230449126182236134239036973933730893336095002F, -0.92387953251128675612818318939678828682241662586364248612F,  -0.55557023301960222474283081394853287437493719075480404593F },
	{ 1,  0.555570233019602224742830813948532874374937190754804045925F, -0.38268343236508977172845998403039886676134456248562704143F,  -0.980785280403230449126182236134239036973933730893336095003F, -0.70710678118654752440084436210484903928483593768847403659F,   0.19509032201612826784828486847702224092769161775195480775F,   0.92387953251128675612818318939678828682241662586364248611F,   0.83146961230254523707878837761790575673856081198724996345F },
	{ 1,  0.195090322016128267848284868477022240927691617751954807755F, -0.923879532511286756128183189396788286822416625863642486115F, -0.555570233019602224742830813948532874374937190754804045925F,  0.707106781186547524400844362104849039284835937688474036587F,  0.831469612302545237078788377617905756738560811987249963448F, -0.382683432365089771728459984030398866761344562485627041431F, -0.980785280403230449126182236134239036973933730893336095004F },
	{ 1, -0.195090322016128267848284868477022240927691617751954807754F, -0.923879532511286756128183189396788286822416625863642486115F,  0.555570233019602224742830813948532874374937190754804045923F,  0.707106781186547524400844362104849039284835937688474036589F, -0.831469612302545237078788377617905756738560811987249963445F, -0.382683432365089771728459984030398866761344562485627041436F,  0.980785280403230449126182236134239036973933730893336095002F },
	{ 1, -0.555570233019602224742830813948532874374937190754804045924F, -0.382683432365089771728459984030398866761344562485627041434F,  0.980785280403230449126182236134239036973933730893336095003F, -0.707106781186547524400844362104849039284835937688474036588F, -0.19509032201612826784828486847702224092769161775195480776F,   0.923879532511286756128183189396788286822416625863642486116F, -0.831469612302545237078788377617905756738560811987249963445F },
	{ 1, -0.831469612302545237078788377617905756738560811987249963446F,  0.382683432365089771728459984030398866761344562485627041434F,  0.195090322016128267848284868477022240927691617751954807755F, -0.707106781186547524400844362104849039284835937688474036589F,  0.980785280403230449126182236134239036973933730893336095003F, -0.923879532511286756128183189396788286822416625863642486115F,  0.55557023301960222474283081394853287437493719075480404592F },
	{ 1, -0.980785280403230449126182236134239036973933730893336095003F,  0.923879532511286756128183189396788286822416625863642486115F, -0.831469612302545237078788377617905756738560811987249963446F,  0.707106781186547524400844362104849039284835937688474036588F, -0.55557023301960222474283081394853287437493719075480404592F,   0.38268343236508977172845998403039886676134456248562704143F,  -0.19509032201612826784828486847702224092769161775195480775F }
};

static const float* matrixA = new float[64] {
	0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273730857504233426880091428756713867187500F,
	0.490392640201615215289621119154617190361022949218750000F,  0.415734806151272617835701339572551660239696502685546875F,  0.277785116509801144335511935423710383474826812744140625F,  0.097545161008064151797469776283833198249340057373046875F, -0.097545161008064096286318545026006177067756652832031250F, -0.277785116509800977802058241650229319930076599121093750F, -0.415734806151272673346852570830378681421279907226562500F, -0.490392640201615215289621119154617190361022949218750000F,
	0.461939766255643369241568052530055865645408630371093750F,  0.191341716182544918645191955874906852841377258300781250F, -0.191341716182544863134040724617079831659793853759765625F, -0.461939766255643369241568052530055865645408630371093750F, -0.461939766255643369241568052530055865645408630371093750F, -0.191341716182545168445372496535128448158502578735351562F,  0.191341716182545001911918802761647384613752365112304688F,  0.461939766255643258219265590014401823282241821289062500F,
	0.415734806151272617835701339572551660239696502685546875F, -0.097545161008064096286318545026006177067756652832031250F, -0.490392640201615215289621119154617190361022949218750000F, -0.277785116509801088824360704165883362293243408203125000F,  0.277785116509800922290907010392402298748493194580078125F,  0.490392640201615215289621119154617190361022949218750000F,  0.097545161008064387719862509129598038271069526672363281F, -0.415734806151272562324550108314724639058113098144531250F,
	0.353553390593273786368655464684707112610340118408203125F, -0.353553390593273730857504233426880091428756713867187500F, -0.353553390593273841879806695942534133791923522949218750F,  0.353553390593273730857504233426880091428756713867187500F,  0.353553390593273841879806695942534133791923522949218750F, -0.353553390593273342279445614622090943157672882080078125F, -0.353553390593273619835201770911226049065589904785156250F,  0.353553390593273286768294383364263921976089477539062500F,
	0.277785116509801144335511935423710383474826812744140625F, -0.490392640201615215289621119154617190361022949218750000F,  0.097545161008064137919681968469376442953944206237792969F,  0.415734806151272784369155033346032723784446716308593750F, -0.415734806151272562324550108314724639058113098144531250F, -0.097545161008064013019591698139265645295381546020507812F,  0.490392640201615326311923581670271232724189758300781250F, -0.277785116509800755757453316618921235203742980957031250F,
	0.191341716182544918645191955874906852841377258300781250F, -0.461939766255643369241568052530055865645408630371093750F,  0.461939766255643258219265590014401823282241821289062500F, -0.191341716182544918645191955874906852841377258300781250F, -0.191341716182545279467674959050782490521669387817382812F,  0.461939766255643369241568052530055865645408630371093750F, -0.461939766255643147196963127498747780919075012207031250F,  0.191341716182544779867313877730339299887418746948242188F,
	0.097545161008064151797469776283833198249340057373046875F, -0.277785116509801088824360704165883362293243408203125000F,  0.415734806151272784369155033346032723784446716308593750F, -0.490392640201615326311923581670271232724189758300781250F,  0.490392640201615215289621119154617190361022949218750000F, -0.415734806151272506813398877056897617876529693603515625F,  0.277785116509800755757453316618921235203742980957031250F, -0.097545161008064276697560046613943995907902717590332031F
};


//  ---------------------------------------------------------------
// |
// |  Naive Sum DCT
// |
//  ---------------------------------------------------------------

void transform8x8_normal(float* input, float* output, size_t width) {
	unsigned char i,j,x,y;
	
	// jump through output
	i = N;
	while (i--) { // outer loop over output
		j = N;
		while (j--) {
			float inner = 0;
			
			// jump through input
			x = N;
			while (x--) { // inner loop over input
				y = N;
				while (y--) {
					inner += input[y + x * width] * cos_2x1iPi_2N[x][i] * cos_2x1iPi_2N[y][j];
				}
			}
			output[j + i * width] = getConstantC(i, j) * inner;
		}
	}
}

void DCT::transform(float* input, float* output, const size_t width, const size_t height) {
	unsigned short x, y;
	const unsigned short width_8 = width / N;
	const size_t blockLineJump = width * N;
	
	y = height / N;
	while (y--) {
		x = width_8;
		while (x--) {
			size_t offset = y * blockLineJump + x * N;
			transform8x8_normal(&input[offset], &output[offset], width);
		}
	}
}



//  ---------------------------------------------------------------
// |
// |  Separated Matrix Multiplication
// |
//  ---------------------------------------------------------------

void multiplyMatrixAWith(float* b, float* result, const size_t width) {
	unsigned char x,y, i;
	y = N;
	while (y--) {
		x = N;
		while (x--) {
			float sum = 0;
			i = N;
			while (i--) { // multiplication loop
				sum += matrixA[y * N + i] * b[i * width + x];
			}
			result[y * N + x] = sum;
		}
	}
}
// two separate loops because its around 300.000 operations per second faster than an if (bool) clause
void multiplyWithTransposedMatrixA(float* a, float* result, const size_t width) {
	unsigned char x,y, i;
	y = N;
	while (y--) {
		x = N;
		while (x--) {
			float sum = 0;
			i = N;
			while (i--) { // multiplication loop
				sum += a[y * N + i] * matrixA[x * N + i];
			}
			result[y * width + x] = sum;
		}
	}
}

void DCT::transform2(float* input, const size_t width, const size_t height) {
	unsigned short x,y;
	const unsigned short width_8 = width / N;
	const size_t blockLineJump = width * N;
	
	float* temp = new float[64];
	
	y = height / N;
	while (y--) {
		x = width_8;
		while (x--) {
			float *ptr = &input[y * blockLineJump + x * N];
			multiplyMatrixAWith(ptr, temp, width); // a * input
			multiplyWithTransposedMatrixA(temp, ptr, width); // input * a^t
		}
	}
	delete[] temp;
}



//  ---------------------------------------------------------------
// |
// |  Invers Calculation
// |
//  ---------------------------------------------------------------

void inverse8x8(float* input, float* output, const size_t width) {
	unsigned char i,j,x,y;
	x = N;
	while (x--) { // outer loop over output
		y = N;
		while (y--) {
			float inner = 0;
			i = N;
			while (i--) { // inner loop over input
				j = N;
				while (j--) {
					float praefix = getConstantC(i, j) * input[j + i * width];
					inner += praefix * cos_2x1iPi_2N[x][i] * cos_2x1iPi_2N[y][j];
				}
			}
			output[y + x * width] = inner;
		}
	}
}

void DCT::inverse(float* input, float* output, const size_t width, const size_t height) {
	unsigned short x, y;
	const unsigned short width_8 = width / N;
	const size_t blockLineJump = width * N;
	
	y = height / N;
	while (y--) {
		x = width_8;
		while (x--) {
			size_t offset = y * blockLineJump + x * N;
			inverse8x8(&input[offset], &output[offset], width);
		}
	}
}
