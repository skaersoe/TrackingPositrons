#ifndef LAUNCHARGUMENTS_H
#define LAUNCHARGUMENTS_H

typedef enum {
  CPU,
  GPU
} launch_device_t;

typedef struct {
  launch_device_t device;
  bool debug;
  int N;
} launch_args_t;

#endif