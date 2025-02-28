import pygame

def gfxdraw_circle(surface, x, y, radius, color):
    """Draw anti-aliased filled circles using pygame.gfxdraw"""
    # Import inside the function to avoid issues if gfxdraw isn't available
    try:
        from pygame import gfxdraw
        
        # Handle both RGB and RGBA colors
        if len(color) == 3:
            r, g, b = color
            a = 255
        else:
            r, g, b, a = color
            
        # Draw filled circle with anti-aliasing
        gfxdraw.filled_circle(surface, x, y, radius, (r, g, b, a))
        gfxdraw.aacircle(surface, x, y, radius, (r, g, b, a))
    except ImportError:
        # Fall back to regular pygame draw if gfxdraw is not available
        pygame.draw.circle(surface, color[:3], (x, y), radius, 0)

def create_grid_pattern(size=20):
    """Create a subtle grid pattern for the background"""
    pattern = pygame.Surface((size, size), pygame.SRCALPHA)
    pattern.fill((0, 0, 0, 0))  # Transparent background
    
    # Draw subtle grid lines
    pygame.draw.line(pattern, (200, 220, 240, 10), (0, 0), (size, 0), 1)
    pygame.draw.line(pattern, (200, 220, 240, 10), (0, 0), (0, size), 1)
    
    return pattern