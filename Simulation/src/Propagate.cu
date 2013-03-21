#include <iostream>
#include <thrust/sort.h>

#include "Simulation/Propagate.h"
#include "Simulation/CudaHelper.h"

namespace Simulation {

  /** Forward declaration */
  __host__ void propagateGPU(simple_particle_t *p, simulator_args_t args);
  __host__ void propagateCPU(simple_particle_t *p, simulator_args_t args);
  __host__ __device__ void timestep(simple_particle_t* p, float dt);
  __global__ void propagateKernel(simple_particle_t *p_, int *keys, kernel_args_t args);

  /** Auxiliary function, intended for debugging. */
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
  void propagate(simple_particle_t *p, simulator_args_t args) {

    args.dt = 0.001;
    args.steps = 2<<12;

    if (args.debug) for (int i=0;i<args.N;i++) print_location(p,i);

    if (args.device == CPU)
      propagateCPU(p,args);
    else
      propagateGPU(p,args);

    if (args.debug) for (int i=0;i<args.N;i++) print_location(p,i);
  }

  /** Regular CPU version for comparison and debugging. */
  __host__
  void propagateCPU(simple_particle_t *p, simulator_args_t args) {
    for (int i=0;i<args.N;i++)
      for (int j=0;j<args.steps;j++)
        timestep(&p[i],args.dt);
  }

  __host__
  void propagateGPU(simple_particle_t *p, simulator_args_t args) {

    // CUDA Parameters
    int threadsPerBlock = 256;
    int blocksPerGrid = (args.N - 1) / threadsPerBlock + 1;
    const unsigned size_particles = args.N*sizeof(simple_particle_t);
    const unsigned size_keys = args.N*sizeof(int);
    const unsigned size_total = size_particles + size_keys;
    kernel_args_t kernel_args = { .N = args.N };

    // Allocate and copy to device
    simple_particle_t *devPtr_p = NULL;
    if (error(cudaMalloc((void**)&devPtr_p,size_total))) return;
    if (error(cudaMemcpy(devPtr_p,p,size_particles,cudaMemcpyHostToDevice))) return;
    int *devPtr_keys = (int*)((long)devPtr_p + (long)size_particles);
    // Thrust wrappers
    thrust::device_ptr<simple_particle_t> devPtr_thrust_p(devPtr_p);
    thrust::device_ptr<int>               devPtr_thrust_keys(devPtr_keys);

    if (args.debug) {
      std::cout << "Copied " << args.N << " instances of size " << sizeof(simple_particle_t) + sizeof(float) << " bytes each, resulting in a total of " << size_total << " bytes of data on the device." << std::endl;
      std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
    }
    // Launch kernel
    for (int i=0; i<10; i++) {
      // Should be dynamic
      kernel_args.steps = 100;
      kernel_args.dt = 0.001;
      propagateKernel<<<blocksPerGrid,threadsPerBlock>>>(devPtr_p,devPtr_keys,kernel_args);
      cudaDeviceSynchronize();
      thrust::sort_by_key(devPtr_thrust_keys, devPtr_thrust_keys + args.N, devPtr_thrust_p);
    }

    // Copy back and free memory
    if (error(cudaMemcpy(p,devPtr_p,size_particles,cudaMemcpyDeviceToHost))) return;
    if (error(cudaFree(devPtr_p))) return;

  }

  /** Called from host function, lets each thread loop over one particle */
  __global__
  void propagateKernel(simple_particle_t *p_, int *keys, kernel_args_t args) {

    int index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index >= args.N) return;

    simple_particle_t* p = &p_[index];

    for (int i=0;i<args.steps;i++) {
      timestep(p,args.dt);
    }

    keys[index] = index % 3;

  }

  /** Applies the magnetic field in a given point */
  __host__ __device__ inline
  void apply_magnetic_field(float *p_src, float *p_tgt, float dt, float q) {
    const float B[] = {10,20,30};
    // Cross product
    p_tgt[0] = dt * (q * (p_src[1] * B[2] - p_src[2] * B[1]) );
    p_tgt[1] = dt * (q * (p_src[2] * B[0] - p_src[0] * B[2]) );
    p_tgt[2] = dt * (q * (p_src[0] * B[1] - p_src[1] * B[0]) );
  }

  /** Auxiliary function for calculating Rung-Kutta coefficients */
  __host__ __device__ inline
  void rk_3d(float *a, float *b, float *c, float d) {
    c[0] = a[0] + d * b[0];
    c[1] = a[1] + d * b[1];
    c[2] = a[2] + d * b[2];
  }

  /** Propagates the particle for one timestep, dt, applying relevant forces */
  __host__ __device__ inline
  void timestep(simple_particle_t* p, float dt) {
    float k1[3], k2[3], k3[3], k4[3], k[3];

    // Rung-Kutta
    apply_magnetic_field(p->p,k1,dt,-1); rk_3d(p->p,k1,k,0.5);
    apply_magnetic_field(k,k2,dt,-1);    rk_3d(p->p,k2,k,0.5);
    apply_magnetic_field(k,k3,dt,-1);    rk_3d(p->p,k3,k,1);
    apply_magnetic_field(k,k4,dt,-1);

    // Update momentum
    p->p[0] += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6;
    p->p[1] += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6;
    p->p[2] += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6;

    // Update position
    p->r[0] += dt * p->p[0];
    p->r[1] += dt * p->p[1];
    p->r[2] += dt * p->p[2];
  }

}