//
//  main.mm
//  MetalTest
//
//  Created by Takahiro Harada on 9/21/16.
//  Copyright © 2016 Takahiro Harada. All rights reserved.
//

#include <iostream>

#import <Metal/Metal.h>

inline
NSString* ec( const char* src )
{
	return [NSString stringWithCString:src encoding:[NSString defaultCStringEncoding]];
}

inline
std::string getKernelString( const char* path )
{
	std::string src;
	FILE* file = fopen(path, "r");
	if( file )
	{
		fseek( file, 0L, SEEK_END );
		size_t binarySize = ftell( file );
		rewind( file );
		char* binary = new char[binarySize];
		fread( binary, sizeof(char), binarySize, file );
		fclose( file );
		src = binary;
		delete [] binary;
	}
	return src;
}

void
saveppm(const char *fname, int w, int h, unsigned char *img)
{
	FILE *fp;
	
	fp = fopen(fname, "wb");
	assert(fp);
	
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", w, h);
	fprintf(fp, "255\n");
	fwrite(img, w * h * 3, 1, fp);
	fclose(fp);
}

#define WIDTH        256
#define HEIGHT       256

int main(int argc, const char * argv[])
{
	id<MTLDevice> device = MTLCreateSystemDefaultDevice();
	
	id<MTLCommandQueue> q = [device newCommandQueue];
	
	const int n = WIDTH*HEIGHT;
	id<MTLBuffer> buff = [device newBufferWithLength:n*sizeof(float)*3 options:MTLResourceOptionCPUCacheModeDefault];
	
	{
		float* d = (float*)[buff contents];
		for(int i=0; i<n*3; i++)
		{
			d[i] = 0.f;
		}
	}
	std::string src = getKernelString("./Kernel.metal");

	NSError *errors;

	id<MTLLibrary> lib = [device newLibraryWithSource:ec(src.c_str()) options:0 error:&errors];

	if( lib == nil )
	{
		NSLog( @"%@",[[errors userInfo] description] );
	}
	id<MTLFunction> func = [lib newFunctionWithName:ec("AoKernel")];
	id<MTLComputePipelineState> state = [device newComputePipelineStateWithFunction:func error:&errors];
	id<MTLCommandBuffer> cmd = [q commandBuffer];
	id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
	[e setComputePipelineState:state];
	
	[e setBuffer:buff offset:0 atIndex:0];
	
	//	run
	const int wgSize = 8;
	MTLSize ng = {WIDTH/wgSize,HEIGHT/wgSize,1};
	MTLSize nt = {wgSize,wgSize,1};
	[e dispatchThreadgroups:ng threadsPerThreadgroup:nt];
	[e endEncoding];
	[cmd commit];
	[cmd waitUntilCompleted];

	{
		unsigned char* dst = new unsigned char[n*3];
		float* d = (float*)[buff contents];
		for(int i=0; i<n*3; i++)
		{
			float x = d[i] * 255.f;
			x = (x >= 255.f)? 255.f:x;
			dst[i] = (unsigned char)x;
		}
		saveppm("ao.ppm", WIDTH, HEIGHT, dst);
		delete [] dst;
	}
	
	std::cout << "Hello Metal Compute!\n";
	
	return 0;
}
