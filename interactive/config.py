import numpy as np

# Display settings
SCREEN_WIDTH = 1280  # Standard HD resolution that fits most screens
SCREEN_HEIGHT = 720  # Standard HD resolution that fits most screens
FPS = 60
BACKGROUND_COLOR = (240, 248, 255)  # Light blue background
SPHERE_COLORS = [(255, 0, 0, 200), (0, 0, 255, 200)]  # Red and blue with alpha

# Sphere settings
SPHERE1_RADIUS = 50  # Radius of the first sphere
SPHERE2_RADIUS = 25  # Radius of the second sphere

# Physics constants
FLUID_VISCOSITY = 100  # Correct original value (was mistakenly changed to 50)
GRAVITY = np.array([0.0, 30.0])  # Gravity direction
TIME_SCALE = .1  # Time scale for more visible motion
MIN_DISTANCE = 2.01  # Minimum allowed distance between sphere centers
WALL_DAMPING = 1.0  # Damping factor for wall collisions
VELOCITY_DAMPING = 1.0  # Reduced damping to allow longer movement
TEMPERATURE = 10000000  # Temperature parameter controlling Brownian motion strength
BROWNIAN_SCALE = 1.0  # Scale factor for Brownian motion forces

# Visual settings
TEXT_BG_COLOR = (30, 30, 30, 180)  # Semi-transparent dark background for text boxes
TEXT_COLOR = (255, 255, 255)  # White text
UI_PADDING = 10  # Padding for UI elements
BOX_CORNER_RADIUS = 8  # Rounded corners for boxes

# Force vector visualization colors
VELOCITY_COLOR = (0, 200, 0)  # Green for velocity vectors
HYDRO_FORCE_COLOR = (255, 165, 0)  # Orange for hydrodynamic forces

# Button styling
BUTTON_BG_COLOR = (60, 60, 70, 220)  # Dark blue-gray with transparency
BUTTON_HOVER_COLOR = (80, 80, 100, 220)  # Lighter when hovered
BUTTON_ACTIVE_COLOR = (100, 180, 100, 220)  # Green when active/enabled
BUTTON_TEXT_COLOR = (255, 255, 255)  # White text
BUTTON_BORDER_COLOR = (120, 120, 140, 255)  # Light border
BUTTON_CORNER_RADIUS = 5  # Rounded corners
