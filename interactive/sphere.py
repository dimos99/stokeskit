import numpy as np
import time
from utils import gfxdraw_circle
import pygame
from config import WALL_DAMPING, TEMPERATURE, BROWNIAN_SCALE

class Sphere:
    def __init__(self, position, radius, mass=1.0, color=(255, 0, 0)):
        self.position = np.array(position, dtype=float)  # Ensure float type
        self.velocity = np.array([0.0, 0.0], dtype=float)
        self.force = np.array([0.0, 0.0], dtype=float)
        self.radius = radius
        self.mass = mass
        self.color = color
        self.dragging = False
        self.start_drag_pos = None
        self.spring_mode = False  # Flag to track if we're in spring mode
        # Visual enhancements
        self.highlight = 0  # For pulsating effect
        self.highlight_direction = 1
        self.last_hydro_force = np.array([0.0, 0.0], dtype=float)  # Store hydrodynamic force for visualization
        
    def apply_force(self, force):
        self.force += force
        
    def update(self, dt):
        # Update position based on velocity
        old_pos = self.position.copy()
        self.position += self.velocity * dt
        
        # Debug output for significant movement
        if np.linalg.norm(self.position - old_pos) > 0.01:
            print(f"Position updated: {old_pos} -> {self.position}, moved: {np.linalg.norm(self.position - old_pos):.3f}")
            
        # Clear forces for next step
        self.force = np.array([0.0, 0.0])
        
    def handle_wall_collision(self, width, height):
        # Left wall
        if self.position[0] - self.radius < 0:
            self.position[0] = self.radius
            self.velocity[0] = -self.velocity[0] * WALL_DAMPING
        
        # Right wall
        elif self.position[0] + self.radius > width:
            self.position[0] = width - self.radius
            self.velocity[0] = -self.velocity[0] * WALL_DAMPING
        
        # Top wall
        if self.position[1] - self.radius < 0:
            self.position[1] = self.radius
            self.velocity[1] = -self.velocity[1] * WALL_DAMPING
            
        # Bottom wall
        elif self.position[1] + self.radius > height:
            self.position[1] = height - self.radius
            self.velocity[1] = -self.velocity[1] * WALL_DAMPING
    
    def draw(self, screen):
        # Enhanced sphere drawing with gradient and lighting effects
        # Main sphere with slight transparency
        base_color = list(self.color[:3])
        
        # Draw the sphere with a gradient effect
        for i in range(int(self.radius), 0, -2):
            # Calculate color gradient from center to edge
            alpha = int(255 * (i / self.radius))
            gradient_color = base_color + [min(alpha, 200)]
            
            # Draw filled circle with antialiasing
            gfxdraw_circle(screen, int(self.position[0]), int(self.position[1]), i, gradient_color)
        
        # Add a highlight effect (specular reflection)
        highlight_radius = int(self.radius * 0.3)
        highlight_pos = (int(self.position[0] - self.radius * 0.3), 
                        int(self.position[1] - self.radius * 0.3))
        # White highlight with transparency
        gfxdraw_circle(screen, highlight_pos[0], highlight_pos[1], highlight_radius, (255, 255, 255, 100))
        
        # Draw outline
        pygame.draw.circle(screen, (base_color[0], base_color[1], base_color[2], 255), 
                          (int(self.position[0]), int(self.position[1])), 
                          int(self.radius), 1)
        
        # Small center point
        pygame.draw.circle(screen, (255, 255, 255), 
                          (int(self.position[0]), int(self.position[1])), 
                          2)
        
    def start_drag(self, pos, spring_mode=False):
        self.dragging = True
        self.start_drag_pos = np.array(pos, dtype=float)  # Ensure float type
        self.velocity = np.array([0.0, 0.0])  # Reset velocity when starting to drag
        self.spring_mode = spring_mode  # Set the spring mode flag
        
    def end_drag(self, pos, dt):
        if self.dragging:
            end_pos = np.array(pos, dtype=float)  # Ensure float type
            # Only calculate velocity if in spring mode
            if self.spring_mode and dt > 0:
                # Scale the throw velocity based on drag distance and direction
                drag_vector = self.start_drag_pos - end_pos
                drag_distance = np.linalg.norm(drag_vector)
                throw_strength = min(drag_distance * 1.0, 800.0)  # Increased strength and max velocity
                self.velocity = drag_vector / dt * 0.2 * throw_strength / drag_distance  # Increased scaling factor
                print(f"Throwing sphere with velocity: {self.velocity}, speed: {np.linalg.norm(self.velocity):.3f}")
            else:
                # In regular mode, just position the sphere with no velocity
                self.velocity = np.array([0.0, 0.0])
            
            self.dragging = False
            self.spring_mode = False
            self.start_drag_pos = None

    def apply_brownian_force(self, dt):
        """Apply random Brownian motion force based on temperature and particle size"""
        # Brownian force variance scales with temperature and inversely with radius
        # Using fluctuation-dissipation theorem: F ~ sqrt(kB*T/dt)
        brownian_amplitude = BROWNIAN_SCALE * TEMPERATURE / (self.radius * np.sqrt(dt))
        
        # Random force with Gaussian distribution in x and y directions
        random_force = np.random.normal(0, brownian_amplitude, 2)
        self.apply_force(random_force)
        
        # Debug output for significant Brownian forces
        if np.linalg.norm(random_force) > 0.5:
            print(f"Brownian force: {random_force}, magnitude: {np.linalg.norm(random_force):.3f}")
        
        return random_force
