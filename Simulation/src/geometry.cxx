#include <cstdio>
#include "TH1.h"


#include <GL/glew.h>

#if defined(WIN32)
#include <GL/wglew.h>
#endif

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include <paramgl.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <helper_cuda.h>
#include <helper_cuda_gl.h>

#include <helper_functions.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <shrQATest.h>

int main (int argc, char const *argv[])
{
    TH1F h("myhitt", "my first histogram!!", 100, 0, 10);
    h.Fill(3232);
    h.Draw("");
    printf("Hello world\n");
    return 0;
}