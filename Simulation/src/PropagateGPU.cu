#include <iostream>
#include <thrust/sort.h>

#include "Simulation/PropagateGPU.cuh"
#include "Simulation/CudaHelper.hh"

namespace na63 {

/** Forward declarations */
__host__ void PropagateGPU(Simulator* simulator);
__host__ __device__ void Step(GPUTrack *t, Float dl);
__global__ void PropagateKernel(KernelPars args);
__host__ cudaError_t DeviceAllocation(Simulator *simulator, KernelPars *p);
__host__ cudaError_t DeviceFree(KernelPars *p);

__host__
void PropagateGPU(Simulator* simulator) {

  if (simulator->debug) {
    for (int i=0;i<5;i++) {
      Track t = simulator->GetTrack(i);
      std::cout << t << std::endl;
    }
  }

  if (simulator->debug) std::cout << "Propagating on GPU..." << std::endl;

  // CUDA Parameters
  int N = simulator->TrackSize();
  int threadsPerBlock = 256;
  int blocksPerGrid = (N - 1) / threadsPerBlock + 1;

  // Initialize arguments
  KernelPars kernel_args;
  kernel_args.N = N;
  kernel_args.tracks = nullptr;
  kernel_args.keys = nullptr;
  kernel_args.materials = nullptr;
  kernel_args.particles = nullptr;
  kernel_args.volume_types = nullptr;
  kernel_args.volumes   = nullptr;

  // Allocate memory on device and copy data
  GPUTrack *tracks = simulator->GPUTracks();
  if (simulator->debug) std::cout << "Generated GPU tracks." << std::endl;
  const unsigned size_tracks = N*sizeof(GPUTrack);
  const unsigned size_keys = N*sizeof(int);
  if (CudaError(cudaMalloc((void**)&kernel_args.tracks,size_tracks))) return;
  if (CudaError(cudaMalloc((void**)&kernel_args.keys,size_keys))) return;
  if (CudaError(cudaMemcpy(kernel_args.tracks,tracks,size_tracks,cudaMemcpyHostToDevice))) return;
  if (simulator->debug) std::cout << "Copied tracks and keys." << std::endl;
  // Thrust wrappers
  thrust::device_ptr<GPUTrack> devptr_thrust_tracks(kernel_args.tracks);
  thrust::device_ptr<int>      devptr_thrust_keys(kernel_args.keys);
  if (CudaError(DeviceAllocation(simulator, &kernel_args))) return;
  if (simulator->debug) std::cout << "Copied geometry." << std::endl;


  if (simulator->debug) {
    std::cout << "Copied " << N << " tracks of size " << sizeof(GPUTrack) << " bytes each, resulting in a total of " << size_tracks << " bytes of data on the device." << std::endl;
    std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
  }

  for (int i=0; i<10; i++) {

    // Should be dynamic
    kernel_args.steps = 100;
    kernel_args.dl = simulator->step_size();

    // Launch kernel
    PropagateKernel<<<blocksPerGrid,threadsPerBlock>>>(kernel_args);
    thrust::sort_by_key(devptr_thrust_keys, devptr_thrust_keys + N, devptr_thrust_tracks);
    cudaDeviceSynchronize();

  } // End kernel launch loop

  // Copy back and free memory
  if (CudaError(cudaMemcpy(tracks,kernel_args.tracks,size_tracks,cudaMemcpyDeviceToHost))) return;
  CudaError(DeviceFree(&kernel_args));

  simulator->CopyBackTracks();

  if (simulator->debug)
    for (int i=0;i<5;i++)
      std::cout << simulator->GetTrack(i) << std::endl;

}

/**
 * Frees used memory on the GPU.
 * @return CUDA error code for first failed instruction, otherwise cudaSuccess.
 */
__host__
cudaError_t DeviceFree(KernelPars *p) {
  cudaError_t err;
  if ((err = cudaFree(p->tracks))       != cudaSuccess) return err;
  if ((err = cudaFree(p->keys))         != cudaSuccess) return err;
  if ((err = cudaFree(p->materials))    != cudaSuccess) return err;
  if ((err = cudaFree(p->particles))    != cudaSuccess) return err;
  if ((err = cudaFree(p->volume_types)) != cudaSuccess) return err;
  if ((err = cudaFree(p->volumes))      != cudaSuccess) return err;
  return cudaSuccess;
}
/**
 * Attempts to allocate room and copy all geometry information to the device.
 * @return CUDA error code for first failed instruction, otherwise cudaSuccess.
 */
__host__
cudaError_t DeviceAllocation(Simulator *simulator, KernelPars *p) {
  cudaError_t err;

  Geometry *geometry = simulator->geometry;

  // Free previously allocated memory
  if (p->materials != nullptr) {
    if ((err = cudaFree(p->materials)) != cudaSuccess) return err;
  }
  if (p->particles != nullptr) {
    if ((err = cudaFree(p->particles)) != cudaSuccess) return err;
  }
  if (p->volume_types != nullptr) {
    if ((err = cudaFree(p->volume_types)) != cudaSuccess) return err;
  }
  if (p->volumes != nullptr) {
    if ((err = cudaFree(p->volumes)) != cudaSuccess) return err;
  }

  // Allocate space
  const unsigned size_material = geometry->materials_size()*sizeof(MaterialPars);
  const unsigned size_particle = simulator->particles_size()*sizeof(ParticlePars);
  const unsigned size_volume_type = geometry->volume_types_size()*sizeof(InsideFunction);
  const unsigned size_volume   = geometry->volumes_size()*sizeof(VolumePars);
  err = cudaMalloc((void**)&p->materials,
                   size_material);
  if (err != cudaSuccess) return err;
  err = cudaMalloc((void**)&p->particles,
                   size_particle);
  if (err != cudaSuccess) return err;
  err = cudaMalloc((void**)&p->volume_types,
                   size_volume_type);
  if (err != cudaSuccess) return err;
  err = cudaMalloc((void**)&p->volumes,
                   size_volume);
  if (err != cudaSuccess) return err;

  // Copy data to device
  err = cudaMemcpy(p->materials,
                   geometry->material_arr(),
                   size_material,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return err;
  err = cudaMemcpy(p->particles,
                   simulator->particle_arr(),
                   size_particle,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return err;
  err = cudaMemcpy(p->volume_types,
                   geometry->volume_type_arr(),
                   size_volume_type,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return err;
  err = cudaMemcpy(p->volumes,
                   geometry->volume_arr(),
                   size_volume,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return err;

  return cudaSuccess;
}

/** Called from host function, lets each thread loop over one particle */
__global__
void PropagateKernel(KernelPars args) {

  int index = threadIdx.x + blockDim.x * blockIdx.x;
  if (index >= args.N) return;

  GPUTrack* t = &args.tracks[index];

  for (int i=0;i<args.steps;i++) {
    Step(t,args.dl);
  }

  args.keys[index] = 0;

}

/** Propagates the particle for one timestep, dt, applying relevant forces */
__host__ __device__ inline
void Step(GPUTrack* t, Float dl) {
  return;
}

} // End namespace na63