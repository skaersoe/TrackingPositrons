#ifndef NA63_SIMULATION_PROCESS_H
#define NA63_SIMULATION_PROCESS_H

namespace na63 {
  
class Process {
public:
  virtual void Query(Track* track, const Material* material, const Float dl) =0;
};

} // End namespace na63

#endif /* NA63_SIMULATION_PROCESS_H */