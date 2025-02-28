import numpy as np
from config import FLUID_VISCOSITY, MIN_DISTANCE

def calculate_hydrodynamic_forces(spheres, selected_sphere=None):
    """Calculate hydrodynamic forces between spheres using StokesKit's resistance matrix."""
    # Skip if any sphere is being dragged
    if selected_sphere is not None:
        return
        
    try:
        # Import stokeskit
        try:
            import stokeskit
        except ImportError:
            print("Error: StokesKit module not found. Using simple drag forces instead.")
            apply_simple_drag_forces(spheres)
            return
            
        # Create arrays for batch computation
        positions = []
        radii = []
        velocities = []
        
        # Collect data from all spheres
        for sphere in spheres:
            # Convert 2D positions to 3D (z=0)
            positions.append(np.array([sphere.position[0], sphere.position[1], 0.0]))
            radii.append(sphere.radius)
            # Convert 2D velocities to 3D (z=0)
            velocities.append(np.array([sphere.velocity[0], sphere.velocity[1], 0.0]))
            
        # Check if spheres are too close - apply minimum spacing
        for i in range(len(spheres)):
            for j in range(i+1, len(spheres)):
                r = positions[j] - positions[i]
                distance = np.linalg.norm(r)
                min_distance = (radii[i] + radii[j]) * MIN_DISTANCE
                
                if distance < min_distance and distance > 1e-10:
                    # Adjust positions to maintain minimum distance
                    r_unit = r / distance
                    adjustment = (min_distance - distance) * 0.5  # Split the adjustment
                    positions[i] -= r_unit * adjustment
                    positions[j] += r_unit * adjustment
        
        # Compute the full resistance matrix (Rinfinity = Minfinity^-1)
        # This accounts for all hydrodynamic interactions
        try:
            # Fix: Use the correct function name from the module
            # Convert Python lists to numpy arrays for the C++ function
            positions_array = np.array(positions, dtype=float)
            radii_array = np.array(radii, dtype=float)
            
            # Try each possible function name that might be available
            try:
                # Option 1: Use computeRinfinityArray from Minfinity submodule
                resistance_matrix = stokeskit.Minfinity.computeRinfinityArray(positions_array, radii_array, FLUID_VISCOSITY)
            except AttributeError:
                try:
                    # Option 2: Use computeRinfinity from main module
                    resistance_matrix = stokeskit.computeRinfinity(positions_array, radii_array, FLUID_VISCOSITY)
                except AttributeError:
                    # Option 3: Use compute_resistance_functions as fallback
                    print("Using scalar resistance functions instead of matrix calculation")
                    # We'll use the simple drag force method instead
                    raise AttributeError("No suitable method found in stokeskit")
                    
        except Exception as e:
            print(f"Resistance matrix calculation failed: {e}")
            apply_simple_drag_forces(spheres)
            return
        
        # Debug info
        print(f"Using full resistance matrix ({resistance_matrix.shape[0]}x{resistance_matrix.shape[1]})")
            
        # Create the correctly sized velocity vector based on matrix shape
        # The matrix size is 11N×11N where N is the number of particles
        # For each particle: 3 for translation, 3 for rotation, 5 for rate of strain
        n_particles = len(spheres)
        total_dof = resistance_matrix.shape[0]  # Should be 11*n_particles
        velocity_vector = np.zeros(total_dof)
        
        # Fill only the translational velocities (first 3 elements per particle)
        for i, v in enumerate(velocities):
            velocity_vector[3*i:3*i+3] = v  # First 3 elements are translational velocities
            # We leave the remaining elements (rotation and strain rate) as zeros
            
        # Calculate forces using resistance matrix: F = R·U
        force_torque_vector = resistance_matrix.dot(velocity_vector)
        
        # Extract only the translational force components for each particle
        for i, sphere in enumerate(spheres):
            # Extract only the forces (first 3 elements per particle)
            force = force_torque_vector[3*i:3*i+3]
            
            # FIX: Invert the force direction to correct the hydrodynamic behavior
            # In Stokes flow, the resistance matrix gives forces that resist motion
            force = -force  # Reverse the sign
            
            # Only use x,y components for 2D simulation
            hydro_force = force[:2]
            sphere.last_hydro_force = hydro_force.copy()  # Store for visualization
            sphere.apply_force(hydro_force)
            
            # Debug forces if significant
            if np.linalg.norm(hydro_force) > 0.1:
                print(f"Sphere {i}: F = {hydro_force}, |F| = {np.linalg.norm(hydro_force):.3f}")
            
    except Exception as e:
        print(f"Hydrodynamic calculation error: {e}")
        import traceback
        traceback.print_exc()
        apply_simple_drag_forces(spheres)

def apply_simple_drag_forces(spheres):
    """Apply simple drag forces as a fallback when StokesKit calculation fails"""
    for sphere in spheres:
        if not sphere.dragging:
            # Simple drag force proportional to velocity in opposite direction
            # Reduced drag coefficient for more visible motion
            drag_coefficient = 6 * np.pi * FLUID_VISCOSITY * sphere.radius * 0.1  # Scaled down Stokes' law
            drag_force = -drag_coefficient * sphere.velocity
            
            # Store the hydrodynamic force for visualization
            sphere.last_hydro_force = drag_force.copy()
            
            sphere.apply_force(drag_force)
            
            # Simple interaction force based on distance
            if len(spheres) >= 2:
                for other in spheres:
                    if sphere is not other:
                        diff = other.position - sphere.position
                        dist = np.linalg.norm(diff)
                        if dist > 0 and dist < 5 * (sphere.radius + other.radius):
                            # Simple repulsion force
                            direction = diff / dist
                            interaction_strength = 2.0 * sphere.radius * other.radius / (dist * dist)  # Increased strength
                            interaction_force = direction * interaction_strength
                            sphere.apply_force(interaction_force)
