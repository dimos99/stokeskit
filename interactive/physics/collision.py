import numpy as np

def handle_sphere_collision(spheres):
    """Handle collision between the spheres"""
    if len(spheres) < 2:
        return
        
    s1, s2 = spheres[0], spheres[1]
    
    # Skip if either sphere is being dragged
    if s1.dragging or s2.dragging:
        return
        
    # Calculate distance between spheres
    diff_vec = s2.position - s1.position
    distance = np.linalg.norm(diff_vec)
    
    # Check for collision
    min_distance = s1.radius + s2.radius
    if distance < min_distance:
        # Move spheres apart
        overlap = min_distance - distance
        direction = diff_vec / distance if distance > 0 else np.array([1.0, 0.0], dtype=float)
        
        # Move proportionally to mass
        total_mass = s1.mass + s2.mass
        s1.position -= direction * overlap * (s2.mass / total_mass)
        s2.position += direction * overlap * (s1.mass / total_mass)
        
        # Elastic collision response
        normal = direction
        relative_velocity = s2.velocity - s1.velocity
        
        # Calculate impulse
        velocity_along_normal = np.dot(relative_velocity, normal)
        if velocity_along_normal > 0:
            return  # Already separating
            
        restitution = 0.8  # Coefficient of restitution
        
        impulse_scalar = -(1 + restitution) * velocity_along_normal
        impulse_scalar /= 1/s1.mass + 1/s2.mass
        
        impulse = impulse_scalar * normal
        
        # Apply impulse
        s1.velocity -= impulse / s1.mass
        s2.velocity += impulse / s2.mass
