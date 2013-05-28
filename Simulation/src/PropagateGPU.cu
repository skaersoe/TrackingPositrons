#include <iostream>
#include <thrust/sort.h>

#include "Simulation/PropagateGPU.cuh"
#include "Simulation/TrackGPU.cuh"
#include "Simulation/CudaHelper.hh"

#include "Simulation/BetheEnergyLossCUDA.cuh"

namespace na63 {

/** Forward declarations */
__host__ void PropagateGPU(Simulator* simulator);
__global__ void PropagateKernel(KernelPars args);
__host__ cudaError_t DeviceAllocation(Simulator *simulator, KernelPars *p);
__host__ cudaError_t DeviceFree(KernelPars *p);

__host__
void PropagateGPU(Simulator* simulator) {

  /*if (simulator->debug) {
    for (int i=0;i<8;i++) {
      Track t = simulator->GetTrack(i);
      std::cout << t << std::endl;
    }
  }*/

  if (simulator->debug) std::cout << "Propagating on GPU..." << std::endl;

  // CUDA Parameters
  int n_tracks_initial = simulator->TrackSize();
  // Use threads equal to the first power of two greater than or equal to the
  // initial amount of tracks.
  int threads = 1 << (int)(log(n_tracks_initial)/log(2));
  if (threads < n_tracks_initial) threads = threads << 1;
  int threadsPerBlock = 256;
  int blocksPerGrid = (threads - 1) / threadsPerBlock + 1;
  // Actually allocate three times as much to allow for pair production
  int N = 3 * n_tracks_initial;

  // Initialize arguments
  KernelPars kernel_args;
  kernel_args.N = N;
  kernel_args.n_volumes = simulator->geometry->volumes_size();
  kernel_args.tracks = nullptr;
  kernel_args.keys = nullptr;
  kernel_args.materials = nullptr;
  kernel_args.particles = nullptr;
  kernel_args.volume_types = nullptr;
  kernel_args.volumes = nullptr;

  // Get track information from simulator
  GPUTrack *tracks = new GPUTrack[N];
  simulator->GPUTracks(tracks);
  // Initialize free tracks
  for (int i=n_tracks_initial;i<N;i++) {
    tracks[i].state = STATE_FREE;
  }
  if (simulator->debug) std::cout << "Generated GPU tracks." << std::endl;

  // Allocate and copy to device
  const unsigned size_tracks_initial = n_tracks_initial * sizeof(GPUTrack);
  const unsigned size_tracks = N*sizeof(GPUTrack);
  const unsigned size_keys = N*sizeof(int);
  if (CudaError(cudaMalloc((void**)&kernel_args.tracks,size_tracks))) return;
  if (CudaError(cudaMalloc((void**)&kernel_args.keys,size_keys))) return;
  if (CudaError(cudaMemcpy(kernel_args.tracks,tracks,size_tracks_initial,cudaMemcpyHostToDevice))) return;
  if (simulator->debug) std::cout << "Copied tracks and keys." << std::endl;
  // Random number generator states
  if (CudaError(cudaMalloc(&kernel_args.rng_states,N*sizeof(curandState)))) return;
  // Thrust wrappers
  thrust::device_ptr<GPUTrack> devptr_thrust_tracks(kernel_args.tracks);
  thrust::device_ptr<int>      devptr_thrust_keys(kernel_args.keys);
  if (CudaError(DeviceAllocation(simulator, &kernel_args))) return;
  if (simulator->debug) std::cout << "Copied geometry." << std::endl;


  if (simulator->debug) {
    std::cout << "Copied " << N << " tracks of size " << sizeof(GPUTrack)
              << " bytes each (of which " << n_tracks_initial 
              << " are initially alive), resulting in a total of " << size_tracks
              << " bytes of data on the device." << std::endl;
    std::cout << "About to initialize " << blocksPerGrid << " blocks of "
              << threadsPerBlock << " threads each, resulting in a total of "
              << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
  }

  int kernel_launches = 0;

  while (1) {

    // Should be dynamic
    kernel_args.steps = 100;
    kernel_args.dl = simulator->step_size;
    kernel_args.rng_seed = kernel_launches;

    // Launch kernel
    kernel_launches++;
    PropagateKernel<<<blocksPerGrid,threadsPerBlock>>>(kernel_args);
    cudaDeviceSynchronize();
    thrust::sort_by_key(devptr_thrust_keys, devptr_thrust_keys + N, devptr_thrust_tracks);
    cudaDeviceSynchronize();

    // Check if propagation is done
    int done = 0;
    if (CudaError(cudaMemcpy((void*)&done,kernel_args.keys,sizeof(int),cudaMemcpyDeviceToHost))) return;
    if (done == TRACK_KEY_DEAD) break;

  } // End kernel launch loop

  if (simulator->debug) std::cout << "Propagation kernel launched "
      << kernel_launches << " times." << std::endl;

  // Copy back and free memory
  if (CudaError(cudaMemcpy(tracks,kernel_args.tracks,size_tracks,cudaMemcpyDeviceToHost))) return;
  CudaError(DeviceFree(&kernel_args));

  simulator->CopyBackTracks(tracks,N);

  /*if (simulator->debug)
    for (int i=0;i<8;i++)
      std::cout << simulator->GetTrack(i) << std::endl;*/

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

  // Set and check thread index
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  if (index >= args.N) return;

  // Get memory pointers
  GPUTrack* track = &args.tracks[index];
  if (track->state != STATE_ALIVE) return;
  ParticlePars* particle = &args.particles[track->particle_index];

  // Initialize random number generator
  curand_init(args.rng_seed,index,0,&args.rng_states[index]);

  // Perform steps
  for (int i=0;i<args.steps;i++) {

    // Update volume if necessary
    track->volume_index = VolumeQuery(*track,args.volumes,args.n_volumes);
    // If out of bounds or out of energy, break
    if (track->volume_index < 0 || track->momentum[3] < particle->mass) {
      track->state = STATE_DEAD;
      break;
    }

    __syncthreads();

    // Propagate position
    Step(*track,*particle,args.dl);

    // Run physics
    if (track->volume_index > 0) {
      const int material_index = args.volumes[track->volume_index].material_index;
      if (material_index >= 0) {
        const MaterialPars *material = &args.materials[args.volumes[track->volume_index].material_index];
        __syncthreads();
        CUDA_BetheEnergyLoss(*track,*particle,*material,args.dl,&args.rng_states[index]);
      }
    }

  }

  __syncthreads();

  // Set key to sort by
  int key = TRACK_KEY_DEAD;
  if (track->state == STATE_ALIVE) {
    key = (int)track->position[0];
  }
  args.keys[index] = key;

}

} // End namespace na63