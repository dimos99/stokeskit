import pygame
import numpy as np
import sys
import time
import math

from config import (SCREEN_WIDTH, SCREEN_HEIGHT, FPS, BACKGROUND_COLOR, SPHERE_COLORS,
                   VELOCITY_DAMPING, TEXT_COLOR, GRAVITY, VELOCITY_COLOR, HYDRO_FORCE_COLOR,
                   SPHERE1_RADIUS, SPHERE2_RADIUS)
from sphere import Sphere
from physics.collision import handle_sphere_collision
from physics.hydrodynamics import calculate_hydrodynamic_forces, apply_simple_drag_forces
from ui.button import Button
from ui.drawing import (draw_vector, draw_legend, draw_text_box, 
                      draw_panel, draw_rounded_rect)
from utils import create_grid_pattern

class HydrodynamicSimulation:
    def __init__(self):
        pygame.init()
        
        # Get actual screen dimensions to ensure window fits
        screen_info = pygame.display.Info()
        avail_width = screen_info.current_w - 100  # Leave margin
        avail_height = screen_info.current_h - 100
        
        # Adjust dimensions if screen is too small
        global SCREEN_WIDTH, SCREEN_HEIGHT
        if SCREEN_WIDTH > avail_width:
            SCREEN_WIDTH = max(800, avail_width)
        if SCREEN_HEIGHT > avail_height:
            SCREEN_HEIGHT = max(600, avail_height)
        
        # Create resizable window to allow user to adjust if needed
        self.screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT), pygame.RESIZABLE)
        pygame.display.set_caption("Hydrodynamic Sphere Interaction")
        self.clock = pygame.time.Clock()
        self.running = True
        
        # Create two spheres using the config values
        self.spheres = [
            Sphere([SCREEN_WIDTH // 4, SCREEN_HEIGHT // 2], SPHERE1_RADIUS, mass=SPHERE1_RADIUS**3, color=SPHERE_COLORS[0]),
            Sphere([3 * SCREEN_WIDTH // 4, SCREEN_HEIGHT // 2], SPHERE2_RADIUS, mass=SPHERE2_RADIUS**3, color=SPHERE_COLORS[1])
        ]
        
        # Font for displaying info - increased sizes for better rendering
        self.font = pygame.font.SysFont('Arial', 16, True)  # Increased size, added antialiasing
        self.title_font = pygame.font.SysFont('Arial', 24, bold=True)  # Increased size
        self.legend_font = pygame.font.SysFont('Arial', 14, True)  # New font for legend
        self.selected_sphere = None
        self.last_update_time = time.time()
        
        # Load background texture
        self.bg_pattern = create_grid_pattern()
        
        # Simulation settings
        self.use_hydrodynamics = True
        self.show_velocity_vectors = True
        self.show_force_vectors = True
        self.use_brownian_motion = True
        self.show_grid = True
        
        # Create buttons - layout in top-left corner
        button_width = 120
        button_height = 30
        button_margin = 10
        button_x = button_margin
        button_y = button_margin
        
        self.buttons = [
            Button(button_x, button_y, button_width, button_height, 
                  "Hydrodynamics", toggle=True, active=self.use_hydrodynamics,
                  tooltip="Toggle between full hydrodynamics and simple drag"),
                  
            Button(button_x, button_y + (button_height + button_margin), button_width, button_height,
                  "Velocity Vectors", toggle=True, active=self.show_velocity_vectors,
                  tooltip="Show/hide velocity vectors"),
                  
            Button(button_x, button_y + 2*(button_height + button_margin), button_width, button_height,
                  "Force Vectors", toggle=True, active=self.show_force_vectors,
                  tooltip="Show/hide hydrodynamic force vectors"),
                  
            Button(button_x, button_y + 3*(button_height + button_margin), button_width, button_height,
                  "Brownian Motion", toggle=True, active=self.use_brownian_motion,
                  tooltip="Enable/disable thermal fluctuations"),
                  
            Button(button_x, button_y + 4*(button_height + button_margin), button_width, button_height,
                  "Show Grid", toggle=True, active=self.show_grid,
                  tooltip="Show/hide background grid"),
                  
            Button(button_x, button_y + 5*(button_height + button_margin), button_width, button_height,
                  "Reset", toggle=False,
                  tooltip="Reset sphere positions and velocities")
        ]
    
    def handle_events(self):
        # Get mouse state
        mouse_pos = pygame.mouse.get_pos()
        mouse_clicked = False
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False
            
            # Handle window resize events
            elif event.type == pygame.VIDEORESIZE:
                global SCREEN_WIDTH, SCREEN_HEIGHT
                SCREEN_WIDTH, SCREEN_HEIGHT = event.size
                self.screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT), pygame.RESIZABLE)
            
            # Track mouse clicks for buttons    
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:  # Left click
                    mouse_clicked = True
                    
                    # Check if clicked on a sphere (only if not clicking on buttons)
                    button_clicked = False
                    for button in self.buttons:
                        if button.rect.collidepoint(mouse_pos):
                            button_clicked = True
                            break
                    
                    if not button_clicked:
                        # Check for sphere clicks
                        spring_mode = pygame.key.get_mods() & pygame.KMOD_SHIFT
                        
                        for i, sphere in enumerate(self.spheres):
                            distance = np.linalg.norm(np.array(mouse_pos) - sphere.position)
                            if distance <= sphere.radius:
                                self.selected_sphere = sphere
                                sphere.start_drag(mouse_pos, spring_mode)
                                break
                
            elif event.type == pygame.MOUSEBUTTONUP:
                if self.selected_sphere:
                    mouse_pos = pygame.mouse.get_pos()
                    dt = time.time() - self.last_update_time
                    self.selected_sphere.end_drag(mouse_pos, dt)
                    self.selected_sphere = None
                    
            elif event.type == pygame.MOUSEMOTION:
                if self.selected_sphere and self.selected_sphere.dragging:
                    self.selected_sphere.position = np.array(pygame.mouse.get_pos(), dtype=float)  # Ensure float type
                
        # Update buttons and handle button clicks
        for i, button in enumerate(self.buttons):
            if button.update(mouse_pos, mouse_clicked):
                # Handle button actions
                if i == 0:  # Hydrodynamics toggle
                    self.use_hydrodynamics = button.active
                elif i == 1:  # Velocity vectors toggle
                    self.show_velocity_vectors = button.active
                elif i == 2:  # Force vectors toggle
                    self.show_force_vectors = button.active
                elif i == 3:  # Brownian motion toggle
                    old_value = self.use_brownian_motion
                    self.use_brownian_motion = button.active
                    # If turning off Brownian motion and user requested immediate stop
                    if old_value and not self.use_brownian_motion and pygame.key.get_mods() & pygame.KMOD_CTRL:
                        self.stop_all_motion()
                elif i == 4:  # Show grid toggle
                    self.show_grid = button.active
                elif i == 5:  # Reset button
                    self.reset_simulation()
    
    def reset_simulation(self):
        """Reset sphere positions and velocities"""
        self.spheres[0].position = np.array([SCREEN_WIDTH // 4, SCREEN_HEIGHT // 2], dtype=float)
        self.spheres[1].position = np.array([3 * SCREEN_WIDTH // 4, SCREEN_HEIGHT // 2], dtype=float)
        self.spheres[0].velocity = np.array([0.0, 0.0], dtype=float)
        self.spheres[1].velocity = np.array([0.0, 0.0], dtype=float)
        self.spheres[0].force = np.array([0.0, 0.0], dtype=float)
        self.spheres[1].force = np.array([0.0, 0.0], dtype=float)
        self.spheres[0].last_hydro_force = np.array([0.0, 0.0], dtype=float)
        self.spheres[1].last_hydro_force = np.array([0.0, 0.0], dtype=float)
    
    def stop_all_motion(self):
        """Stop all motion of spheres immediately"""
        for sphere in self.spheres:
            sphere.velocity = np.array([0.0, 0.0], dtype=float)
            sphere.force = np.array([0.0, 0.0], dtype=float)
            sphere.last_hydro_force = np.array([0.0, 0.0], dtype=float)

    def update(self):
        current_time = time.time()
        dt = (current_time - self.last_update_time) * VELOCITY_DAMPING
        self.last_update_time = current_time
        
        # Skip if dt is too large (e.g., after dragging)
        if dt > 0.1:
            dt = 0.1  # Increased from 0.01 to 0.1
        
        try:
            # Calculate hydrodynamic forces based on setting
            if self.use_hydrodynamics:
                calculate_hydrodynamic_forces(self.spheres, self.selected_sphere)
            else:
                apply_simple_drag_forces(self.spheres)
            
            # Apply gravity and Brownian motion to both spheres
            for sphere in self.spheres:
                if not sphere.dragging:
                    # Apply gravity
                    gravity_force = GRAVITY * sphere.mass
                    sphere.apply_force(gravity_force)
                    
                    # Apply Brownian motion forces if enabled
                    if self.use_brownian_motion:
                        brownian_force = sphere.apply_brownian_force(dt)
                    
                    # Show total forces being applied
                    if np.linalg.norm(sphere.force) > 0.1:
                        print(f"Applying force: {sphere.force}, magnitude: {np.linalg.norm(sphere.force):.3f}")
                    
                    # Update velocity based on force
                    old_vel = sphere.velocity.copy()
                    sphere.velocity += sphere.force / sphere.mass * dt
                    
                    # Show velocity changes
                    if np.linalg.norm(sphere.velocity - old_vel) > 0.01:
                        print(f"Velocity updated: {old_vel} -> {sphere.velocity}, delta: {np.linalg.norm(sphere.velocity - old_vel):.3f}")
                    
                    # Apply damping (fluid drag) - reduced for more visible motion
                    sphere.velocity *= VELOCITY_DAMPING
                    
                    # Update position
                    sphere.update(dt)
                    
                    # Debug output for velocity
                    if np.linalg.norm(sphere.velocity) > 0.1:
                        print(f"Sphere moving with velocity: {sphere.velocity}, speed: {np.linalg.norm(sphere.velocity):.3f}")
                    
                    # Handle wall collisions
                    sphere.handle_wall_collision(SCREEN_WIDTH, SCREEN_HEIGHT)
            
            # Handle collisions between spheres
            handle_sphere_collision(self.spheres)
        except Exception as e:
            print(f"Error in simulation update: {e}")
            import traceback
            traceback.print_exc()
    
    def draw(self):
        # Fill with background color
        self.screen.fill(BACKGROUND_COLOR)
        
        # Draw background grid pattern if enabled
        if self.show_grid:
            for x in range(0, SCREEN_WIDTH, 20):
                for y in range(0, SCREEN_HEIGHT, 20):
                    self.screen.blit(self.bg_pattern, (x, y))
                
        # Draw a subtle water-like effect at the bottom
        water_height = int(SCREEN_HEIGHT * 0.1)  # 10% of screen height
        water_surface = pygame.Surface((SCREEN_WIDTH, water_height), pygame.SRCALPHA)
        
        # Create a gradient from transparent to light blue
        for y in range(water_height):
            alpha = int(50 * y / water_height)  # Gradually increase alpha
            pygame.draw.line(water_surface, (100, 150, 255, alpha), 
                           (0, y), (SCREEN_WIDTH, y))
        
        self.screen.blit(water_surface, (0, SCREEN_HEIGHT - water_height))
        
        # Draw spheres
        for sphere in self.spheres:
            sphere.draw(self.screen)
            
            if not sphere.dragging:
                # Draw velocity vector if enabled
                if self.show_velocity_vectors:
                    velocity_magnitude = np.linalg.norm(sphere.velocity)
                    draw_vector(self.screen, sphere.position, sphere.velocity, VELOCITY_COLOR, 5.0, "V")
                    
                    # Display velocity magnitude
                    if velocity_magnitude > 0.001:
                        font = pygame.font.SysFont('Arial', 10)
                        text = font.render(f"v={velocity_magnitude:.2f}", True, VELOCITY_COLOR)
                        text_pos = (int(sphere.position[0]) - 40, int(sphere.position[1]) - 15)
                        self.screen.blit(text, text_pos)
                
                # Draw hydrodynamic force vector if enabled
                if self.show_force_vectors:
                    # Linear scaling for forces that preserves the proportional relationship with velocity
                    # In Stokes flow: F = 6πηRv (force is proportional to velocity)
                    force_magnitude = np.linalg.norm(sphere.last_hydro_force)
                    
                    # Only draw if there's a meaningful force
                    if force_magnitude > 0.001:
                        # Use direct linear scaling to maintain proportionality with velocity
                        draw_vector(self.screen, sphere.position, sphere.last_hydro_force, HYDRO_FORCE_COLOR, 1E-3, "F")
                        
                        # Display force magnitude near the vector
                        font = pygame.font.SysFont('Arial', 10)
                        text = font.render(f"F={force_magnitude:.2f}", True, HYDRO_FORCE_COLOR)
                        text_pos = (int(sphere.position[0]) + 25, int(sphere.position[1]) + 15)
                        self.screen.blit(text, text_pos)
            
            # Draw spring indicator with enhanced visuals
            elif sphere.spring_mode and sphere.start_drag_pos is not None:
                # Draw a spring-like line
                start_pos = tuple(sphere.position)
                end_pos = tuple(sphere.start_drag_pos)
                
                # Calculate spring segments
                dist = np.linalg.norm(sphere.position - sphere.start_drag_pos)
                segments = max(5, min(int(dist / 10), 20))
                direction = (sphere.start_drag_pos - sphere.position) / dist
                perp = np.array([-direction[1], direction[0]])
                amplitude = min(dist * 0.1, 10)
                
                # Draw spring segments
                prev_point = start_pos
                for i in range(1, segments+1):
                    t = i / segments
                    pos = sphere.position + direction * dist * t
                    
                    # Add sine wave perpendicular offset for spring effect
                    if i > 0 and i < segments:
                        offset = perp * amplitude * math.sin(i * math.pi)
                        pos = pos + offset
                    
                    point = tuple(pos)
                    pygame.draw.line(self.screen, (255, 50, 50), prev_point, point, 2)
                    prev_point = point
        
        # Draw title in a clean box at the top
        title = "Hydrodynamic Sphere Interaction"
        title_y = 30
        draw_text_box(self.screen, title, SCREEN_WIDTH // 2, title_y, self.title_font, centered=True)
        
        # Draw legend if vectors are displayed
        if self.show_velocity_vectors or self.show_force_vectors:
            legend_items = []
            if self.show_velocity_vectors:
                legend_items.append(("Velocity (V)", VELOCITY_COLOR))
            if self.show_force_vectors:
                legend_items.append(("Hydro Force (F)", HYDRO_FORCE_COLOR))
                
            if legend_items:
                legend_x = SCREEN_WIDTH - 200
                legend_y = 20
                draw_legend(self.screen, legend_x, legend_y, legend_items, self.legend_font)
        
        # Draw info panel with fixed dimensions that fit the text
        info_text = [
            "Drag: move sphere",
            "Shift+Drag: spring action",
            "Ctrl+Click Brownian toggle: stop motion",
            f"Sphere 1: pos=({self.spheres[0].position[0]:.0f}, {self.spheres[0].position[1]:.0f})",
            f"Sphere 2: pos=({self.spheres[1].position[0]:.0f}, {self.spheres[1].position[1]:.0f})"
        ]
        
        # Increased panel size for better text rendering
        panel_width = 300
        panel_height = len(info_text) * 22 + 20  # Increased height for text
        panel_x = 15
        panel_y = SCREEN_HEIGHT - panel_height - 15
        
        draw_panel(self.screen, panel_x, panel_y, panel_width, panel_height)
        
        # Draw text with more padding and antialiasing
        for i, text in enumerate(info_text):
            text_surface = self.font.render(text, True, TEXT_COLOR)
            self.screen.blit(text_surface, (panel_x + 10, panel_y + 10 + i * 22))
            
        # Draw buttons
        for button in self.buttons:
            button.draw(self.screen)
        
        pygame.display.flip()
    
    def run(self):
        while self.running:
            self.handle_events()
            self.update()
            self.draw()
            self.clock.tick(FPS)
        
        pygame.quit()
        sys.exit()
