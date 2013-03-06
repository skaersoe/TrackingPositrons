#include <iostream>

#include "Simulation/Propagate.h"
#include "Simulation/CudaHelper.h"

__host__ __device__ void timestep(simple_particle_t* p, float dt);
__host__   void propagateCPU(simple_particle_t *p, const int N, const int steps, float dt);
__global__ void propagateGPU(simple_particle_t *p, const int N, const int steps, float dt);

/**
 * Auxiliary function, intended for debugging.
 */
__host__
void print_location(simple_particle_t *p, const int i) {
  std::cout << "Particle " << i << " is at (" << p[i].r[0] << "," << p[i].r[1]
    << "," << p[i].r[2] << ") with momentum (" << p[i].p[0] << "," << p[i].p[1]
    << "," << p[i].p[2] << ")" << std::endl;
}

/**
 * Controller function. Propagates the particles on the CPU or the GPU depending
 * on the input arguments.
 */
__host__
void propagate(simple_particle_t *p, const launch_args_t args) {

  /* Global constants */
  const int N = args.N;
  float dt = 1e-11; // 10ps
  const int steps = (int)((10./2.99792e8)/dt); // (100m / c) / dt

  if (args.debug) print_location(p,5);

  if (args.device == CPU) {

    propagateCPU(p,N,steps,dt);

  } else

  if (args.device == GPU) {

    /* CUDA Parameters */
    int threadsPerBlock = 256;
    int blocksPerGrid = (N - 1) / threadsPerBlock + 1;
    const int dataSize = N*sizeof(simple_particle_t);

    /* Allocate and copy to device */
    simple_particle_t* devicePtr = NULL;
    
    if (error(cudaMalloc((void**)&devicePtr,dataSize))) return;
    if (error(cudaMemcpy((void*)devicePtr,p,dataSize,cudaMemcpyHostToDevice))) return;

    if (args.debug) {
      std::cout << "Copied " << N << " instances of size " << sizeof(simple_particle_t) << " bytes each, resulting in a total of " << dataSize << " bytes of data on the device." << std::endl;
      std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
    }
    /* Launch kernel */
    propagateGPU<<<blocksPerGrid,threadsPerBlock>>>(devicePtr,N,steps,dt);
    cudaDeviceSynchronize();

    /* Copy back and free memory */
    if (error(cudaMemcpy(p,devicePtr,dataSize,cudaMemcpyDeviceToHost))) return;
    if (error(cudaFree(devicePtr))) return;

  } // End GPU Propagataion

  if (args.debug) print_location(p,5);
}

/**
 * Regular CPU version for comparison and debugging.
 */
__host__
void propagateCPU(simple_particle_t *p, const int N, const int steps, float dt) {
  for (int i=0;i<N;i++)
    for (int j=0;j<steps;j++)
      timestep(&p[i],dt);
}

/**
 * Called from host function, lets each thread loop over one particle
 */
__global__
void propagateGPU(simple_particle_t *p_, const int N, const int steps, float dt) {
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  if (index >= N) return;
  simple_particle_t* p = &p_[index];

  for (int j=0;j<steps;j++)
    timestep(p,dt);
}

/**
 * Applies the magnetic field in a given point
 */
__host__ __device__
void apply_magnetic_field(float *p_src, float *p_tgt, float dt, float q) {
  const float B[] = {0,0,666};
  // Cross product
  p_tgt[0] = dt * (q * (p_src[1] * B[2] - p_src[2] * B[1]) );
  p_tgt[1] = dt * (q * (p_src[2] * B[0] - p_src[0] * B[2]) );
  p_tgt[2] = dt * (q * (p_src[0] * B[1] - p_src[1] * B[0]) );
}

/**
 * Auxiliary function for calculating Rung-Kutta coefficients
 */
__host__ __device__
void rk_3d(float *a, float *b, float *c, float d) {
  c[0] = a[0] + d * b[0];
  c[1] = a[1] + d * b[1];
  c[2] = a[2] + d * b[2];
}

/**
 * Propagates the particle for one timestep, dt, applying relevant forces
 */
__host__ __device__
void timestep(simple_particle_t* p, float dt) {
  float k1[3], k2[3], k3[3], k4[3], k[3];

  // Rung-Kutta
  apply_magnetic_field(p->p,k1,dt,p->q); rk_3d(p->p,k1,k,0.5);
  apply_magnetic_field(k,k2,dt,p->q);    rk_3d(p->p,k2,k,0.5);
  apply_magnetic_field(k,k3,dt,p->q);    rk_3d(p->p,k3,k,1);
  apply_magnetic_field(k,k4,dt,p->q);

  // Update momentum
  p->p[0] += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6;
  p->p[1] += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6;
  p->p[2] += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6;

  // Update position
  p->r[0] += dt * p->p[0];
  p->r[1] += dt * p->p[1];
  p->r[2] += dt * p->p[2];
}