# Interactive Hydrodynamic Spheres Simulation

This project demonstrates hydrodynamic interactions between two spheres in a viscous fluid using StokesKit for the underlying physics calculations.

## Description

This interactive environment allows you to manipulate two spheres and observe their hydrodynamic interactions in real-time. The simulation uses:

- **StokesKit**: For hydrodynamic calculations (resistance functions)
- **Pygame**: For visualization and user interaction

## Features

- Interactive dragging and throwing of spheres
- Spring-like behavior with Shift+Drag
- Real-time visualization of hydrodynamic interactions
- Toggleable display of velocity and force vectors
- Brownian motion simulation
- Physically accurate fluid dynamics based on Stokes flow
- Elastic collisions between spheres and boundaries
- Adjustable parameters through UI buttons
- Visual grid background with toggle option
- Antialiased sphere rendering with gradient and lighting effects

## Requirements

- Python 3.6+
- Pygame 2.0+
- NumPy 1.19+
- StokesKit

## Installation

1. Make sure StokesKit is properly installed and accessible in your Python environment
2. Install additional required packages:

```bash
pip install pygame numpy
```

## How to Run

Run the interactive simulation using the main script:

```bash
python main.py
```

## Project Structure

The project is organized as follows:

- `/interactive/` - Contains all the interactive simulation code
  - `main.py` - Entry point script
  - `simulation.py` - Main simulation class
  - `sphere.py` - Sphere class implementation
  - `config.py` - Configuration settings
  - `utils.py` - Utility functions
  - `/physics/` - Physics calculations
    - `collision.py` - Collision handling between spheres
    - `hydrodynamics.py` - Hydrodynamic forces calculation
  - `/ui/` - User interface components
    - `button.py` - Interactive button implementation
    - `drawing.py` - Drawing utilities for vectors, panels, etc.

## Controls

### Basic Controls:
- **Click and drag** a sphere to move it
- **Release** the mouse button to place the sphere with no velocity
- **Shift+Click and drag** to create a spring-like effect and throw the sphere
- **Ctrl+Click** the Brownian Motion button to stop all motion immediately

### UI Controls:
- **Hydrodynamics** - Toggle between full hydrodynamics and simple drag
- **Velocity Vectors** - Show/hide velocity vectors
- **Force Vectors** - Show/hide hydrodynamic force vectors
- **Brownian Motion** - Enable/disable thermal fluctuations
- **Show Grid** - Toggle background grid visibility
- **Reset** - Reset sphere positions and velocities

## Physics Details

The simulation uses StokesKit to calculate hydrodynamic interactions between spheres in a Stokes flow regime:

1. The resistance functions are used to calculate forces between spheres
2. Forces are applied to spheres based on their relative velocities and positions
3. Sphere motion includes gravity, hydrodynamic forces, Brownian motion, and collision responses
4. The simulation demonstrates key concepts in low Reynolds number fluid dynamics:
   - Long-range hydrodynamic interactions
   - Viscous drag effects
   - Thermal fluctuations (Brownian motion)

## Advanced Features

### Vector Visualization
The simulation provides real-time visualization of:
- Velocity vectors (green)
- Hydrodynamic force vectors (orange)

### Physical Accuracy
- Properly scaled hydrodynamic interactions using resistance matrices
- Mass-proportional gravity effects
- Viscous damping of motion
- Elastic collisions with conservation of momentum
- Brownian motion scaled by temperature and particle size

## Troubleshooting

If you encounter performance issues:
- Disable Brownian motion (using the button)
- Reduce the display resolution in config.py
- Simplify the physics model by turning off hydrodynamics
- Ensure your Python environment has access to StokesKit
