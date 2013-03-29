#include <iostream>
#include <thrust/sort.h>

#include "Simulation/Propagate.cuh"
#include "Simulation/CudaHelper.hh"

namespace na63 {

  /** Forward declaration */
  __host__ void PropagateGPU(Track *t, SimulatorPars args);
  __host__ void PropagateCPU(Track *t, SimulatorPars args);
  __host__ __device__ void Timestep(Track *t, float dt);
  __global__ void PropagateKernel(Track *t_, int *keys, KernelPars args);
  __host__ cudaError_t AllocateGeometry(Geometry *geometry, KernelPars *p);
  __host__ cudaError_t FreeGeometry(KernelPars *p);

  /** Auxiliary function, intended for debugging. */
  __host__
  void PrintLocation(Track *t, const int i) {
    std::cout << "Particle " << i << " is at (" << t[i].r[0] << "," << t[i].r[1]
      << "," << t[i].r[2] << ") with momentum (" << t[i].p[0] << "," << t[i].p[1]
      << "," << t[i].p[2] << ")" << std::endl;
  }

  /**
   * Controller function. Propagates the particles on the CPU or the GPU depending
   * on the input arguments.
   */
   __host__
  void Propagate(Track *t, SimulatorPars args) {

    if (args.debug) for (int i=0;i<5;i++) PrintLocation(t,i);

    if (args.device == CPU)
      PropagateCPU(t,args);
    else
      PropagateGPU(t,args);

    if (args.debug) for (int i=0;i<5;i++) PrintLocation(t,i);
  }

  /** Regular CPU version for comparison and debugging. */
  __host__
  void PropagateCPU(Track *t, SimulatorPars args) {
    // Not currently implemented
    t = t;
    args = args;
  }

  __host__
  void PropagateGPU(Track *tracks, SimulatorPars args) {

    // CUDA Parameters
    int threadsPerBlock = 256;
    int blocksPerGrid = (args.N - 1) / threadsPerBlock + 1;
    const unsigned size_tracks = args.N*sizeof(Track);
    const unsigned size_keys = args.N*sizeof(int);
    const unsigned size_total = size_tracks + size_keys;
    KernelPars kernel_args;
    kernel_args.N = args.N;
    kernel_args.material_arr = NULL;
    kernel_args.particle_arr = NULL;
    kernel_args.volume_arr   = NULL;

    // Allocate geometry arrays (might be dynamic later)
    if (CudaError(AllocateGeometry(args.geometry, &kernel_args))) return;

    // Allocate tracks and copy to device
    Track *devptr_tracks = NULL;
    if (CudaError(cudaMalloc((void**)&devptr_tracks,size_total))) return;
    if (CudaError(cudaMemcpy(devptr_tracks,tracks,size_tracks,cudaMemcpyHostToDevice))) return;
    int *devptr_keys = (int*)((long)devptr_tracks + (long)size_tracks);
    // Thrust wrappers
    thrust::device_ptr<Track> devptr_thrust_tracks(devptr_tracks);
    thrust::device_ptr<int>   devptr_thrust_keys(devptr_keys);

    if (args.debug) {
      std::cout << "Copied " << args.N << " instances of size " << sizeof(Track) + sizeof(float) << " bytes each, resulting in a total of " << size_total << " bytes of data on the device." << std::endl;
      std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
    }
    // Launch kernel
    for (int i=0; i<10; i++) {
      // Should be dynamic //
      kernel_args.steps = 100;
      kernel_args.dt = 0.001;
      // ----------------- //
      PropagateKernel<<<blocksPerGrid,threadsPerBlock>>>(devptr_tracks,devptr_keys,kernel_args);
      cudaDeviceSynchronize();
      thrust::sort_by_key(devptr_thrust_keys, devptr_thrust_keys + args.N, devptr_thrust_tracks);
    } // End kernel launch loop

    // Copy back and free memory
    if (CudaError(cudaMemcpy(tracks,devptr_tracks,size_tracks,cudaMemcpyDeviceToHost))) return;
    if (CudaError(cudaFree(devptr_tracks))) return;
    if (CudaError(FreeGeometry(&kernel_args))) return;

  }

  __host__
  cudaError_t FreeGeometry(KernelPars *p) {
    cudaError_t err;
    if ((err = cudaFree(p->material_arr)) != cudaSuccess) return err;
    if ((err = cudaFree(p->particle_arr)) != cudaSuccess) return err;
    if ((err = cudaFree(p->volume_arr))   != cudaSuccess) return err;
    return cudaSuccess;
  }
  __host__
  cudaError_t AllocateGeometry(Geometry *geometry, KernelPars *p) {
    cudaError_t err;

    // Free previously allocated geometry if necessary
    if (p->material_arr != NULL) {
      if ((err = cudaFree(p->material_arr)) != cudaSuccess) return err;
    }
    if (p->particle_arr != NULL) {
      if ((err = cudaFree(p->particle_arr)) != cudaSuccess) return err;
    }
    if (p->volume_arr != NULL) {
      if ((err = cudaFree(p->volume_arr)) != cudaSuccess) return err;
    }

    // Allocate space
    const unsigned size_material = geometry->materials_size()*sizeof(MaterialPars);
    const unsigned size_particle = geometry->particles_size()*sizeof(ParticlePars);
    const unsigned size_volume   = geometry->volumes_size()*sizeof(VolumePars);
    err = cudaMalloc((void**)&p->material_arr,
                     size_material);
    if (err != cudaSuccess) return err;
    err = cudaMalloc((void**)&p->particle_arr,
                     size_particle);
    if (err != cudaSuccess) return err;
    err = cudaMalloc((void**)&p->volume_arr,
                     size_volume);
    if (err != cudaSuccess) return err;

    // Copy geometry to device (could potentially be asynchronous?)
    err = cudaMemcpy(p->material_arr,
                     geometry->material_arr(),
                     size_material,
                     cudaMemcpyHostToDevice);
    if (err != cudaSuccess) return err;
    err = cudaMemcpy(p->particle_arr,
                     geometry->particle_arr(),
                     size_particle,
                     cudaMemcpyHostToDevice);
    if (err != cudaSuccess) return err;
    err = cudaMemcpy(p->volume_arr,
                     geometry->volume_arr(),
                     size_volume,
                     cudaMemcpyHostToDevice);
    if (err != cudaSuccess) return err;

    return cudaSuccess;
  }

  /** Called from host function, lets each thread loop over one particle */
  __global__
  void PropagateKernel(Track *t_, int *keys, KernelPars args) {

    int index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index >= args.N) return;

    Track* t = &t_[index];

    for (int i=0;i<args.steps;i++) {
      Timestep(t,args.dt);
    }

    keys[index] = index % 3;

  }

  /** Applies the magnetic field in a given point */
  __host__ __device__ inline
  void ApplyMagneticField(float *t_src, float *t_tgt, float dt, float q) {
    const float B[] = {10,20,30};
    // Cross product
    t_tgt[0] = dt * (q * (t_src[1] * B[2] - t_src[2] * B[1]) );
    t_tgt[1] = dt * (q * (t_src[2] * B[0] - t_src[0] * B[2]) );
    t_tgt[2] = dt * (q * (t_src[0] * B[1] - t_src[1] * B[0]) );
  }

  /** Auxiliary function for calculating Rung-Kutta coefficients */
  __host__ __device__ inline
  void RK3D(float *a, float *b, float *c, float d) {
    c[0] = a[0] + d * b[0];
    c[1] = a[1] + d * b[1];
    c[2] = a[2] + d * b[2];
  }

  /** Propagates the particle for one timestep, dt, applying relevant forces */
  __host__ __device__ inline
  void Timestep(Track* t, float dt) {
    float k1[3], k2[3], k3[3], k4[3], k[3];

    // Rung-Kutta
    ApplyMagneticField(t->p,k1,dt,-1); RK3D(t->p,k1,k,0.5);
    ApplyMagneticField(k,k2,dt,-1);    RK3D(t->p,k2,k,0.5);
    ApplyMagneticField(k,k3,dt,-1);    RK3D(t->p,k3,k,1);
    ApplyMagneticField(k,k4,dt,-1);

    // Update momentum
    t->p[0] += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6;
    t->p[1] += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6;
    t->p[2] += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6;

    // Update position
    t->r[0] += dt * t->p[0];
    t->r[1] += dt * t->p[1];
    t->r[2] += dt * t->p[2];
  }

}