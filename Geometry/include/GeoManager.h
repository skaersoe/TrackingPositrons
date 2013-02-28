#ifndef GEOMANAGER_H_VO05Y572
#define GEOMANAGER_H_VO05Y572

/*
    A GeoManager instance handles the user access to the geometry model and 
    is the point of registration of objects constituting the geometry.
    
    This method should give the user access to the following:
    
    - Setup:
        * Defining the overall geometrical frame
        * Adding detector and material objects
        * Adding magnetic fields in the form of LUTs or parameterisations
        
    - Query:
        * position dependent lookup of material
        * position dependent lookup of magnetic fields
        * position dependent lookup of detectors
        * Calculate radiation length between two positions
        * 
    
    - Geometry services
        * Allow query in arbitrary coordinate systems
        * Place objects in rotated, translated and mirroed states
        * Allow for ROOT EvE and Geometry objects to be returned for 
          event viewers.
          
*/
class GeoManager {
public:
    GeoManager();
    
    
protected:
private:
};

#endif /* end of include guard: GEOMANAGER_H_VO05Y572 */
