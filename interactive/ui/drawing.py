import pygame
import numpy as np
from config import BOX_CORNER_RADIUS, TEXT_BG_COLOR, TEXT_COLOR, UI_PADDING

def draw_vector(screen, origin, vector, color, scale=1.0, label=None):
    """Draw a vector arrow from origin in the direction of vector"""
    if np.linalg.norm(vector) < 0.05:
        return  # Don't draw very small vectors
    
    # Scale and normalize
    magnitude = np.linalg.norm(vector)
    length = min(magnitude * scale, 150)  # Increased max length
    
    if magnitude > 0:
        direction = vector / magnitude
        arrow_perp = np.array([-direction[1], direction[0]])
    else:
        return
    
    # Calculate end point
    end_pos = origin + direction * length
    
    # Draw line with thickness proportional to magnitude
    thickness = max(2, min(int(magnitude/5), 6))  # Adjusted thickness scaling
    pygame.draw.line(screen, color, 
                   (int(origin[0]), int(origin[1])),
                   (int(end_pos[0]), int(end_pos[1])), 
                   thickness)
    
    # Draw arrow head if vector is long enough
    if length > 8:  # Reduced minimum length for arrow head
        # Calculate arrow head points - size proportional to vector length
        arrow_size = min(8 + length/8, 25)  # Adjusted scaling
        
        arrow_p1 = end_pos - direction * arrow_size + arrow_perp * arrow_size * 0.4
        arrow_p2 = end_pos - direction * arrow_size - arrow_perp * arrow_size * 0.4
        
        pygame.draw.polygon(screen, color, 
                          [tuple(map(int, end_pos)), 
                           tuple(map(int, arrow_p1)), 
                           tuple(map(int, arrow_p2))])
    
    # Add label if provided - improved positioning
    if label:
        label_offset = 15  # Increased offset
        mid_point = origin + direction * (length * 0.5)
        offset = arrow_perp * label_offset
        label_pos = mid_point + offset
        
        font = pygame.font.SysFont('Arial', 14, bold=True)  # Increased size
        text = font.render(label, True, color)
        text_rect = text.get_rect(center=tuple(map(int, label_pos)))
        screen.blit(text, text_rect)

def draw_legend(screen, x, y, items, legend_font):
    """Draw a legend with vector symbols and labels"""
    line_length = 25  # Increased length
    spacing = 30  # Increased spacing
    padding = 10  # Increased padding
    
    # Calculate total height and width
    total_height = len(items) * spacing + padding * 2
    max_text_width = max(legend_font.render(item[0], True, TEXT_COLOR).get_width() 
                       for item in items)
    total_width = line_length + 15 + max_text_width + padding * 2
    
    # Draw legend background
    legend_rect = pygame.Rect(x, y, total_width, total_height)
    draw_rounded_rect(screen, legend_rect, TEXT_BG_COLOR)
    
    # Draw each item
    for i, (label, color) in enumerate(items):
        # Draw line
        line_y = y + padding + i * spacing + spacing // 2
        pygame.draw.line(screen, color, 
                       (x + padding, line_y), 
                       (x + padding + line_length, line_y), 
                       2)
        
        # Draw arrowhead
        pygame.draw.polygon(screen, color, 
                         [(x + padding + line_length, line_y),
                          (x + padding + line_length - 6, line_y - 4),
                          (x + padding + line_length - 6, line_y + 4)])
        
        # Draw label with improved font
        text = legend_font.render(label, True, TEXT_COLOR)
        screen.blit(text, (x + padding + line_length + 8, line_y - 8))

def draw_text_box(screen, text, x, y, font, centered=False, width=None):
    """Draw text in a semi-transparent box"""
    # Render text to get its dimensions
    text_surface = font.render(text, True, TEXT_COLOR)
    text_rect = text_surface.get_rect()
    
    if centered:
        text_rect.center = (x, y)
    else:
        text_rect.topleft = (x, y)
        
    # Create box rect with conservative padding
    box_rect = text_rect.inflate(UI_PADDING * 3, UI_PADDING * 1.5)
    
    if width and width > box_rect.width:
        box_rect.width = width
        
    # Ensure box is properly positioned
    if centered:
        box_rect.center = (x, y)
        text_rect.center = box_rect.center
    
    # Draw rounded box background
    draw_rounded_rect(screen, box_rect, TEXT_BG_COLOR)
        
    # Draw text
    screen.blit(text_surface, text_rect)
    
    return box_rect
    
def draw_panel(screen, x, y, width, height):
    """Draw a semi-transparent panel with rounded corners"""
    panel_rect = pygame.Rect(x, y, width, height)
    draw_rounded_rect(screen, panel_rect, TEXT_BG_COLOR)
    
def draw_rounded_rect(screen, rect, color, border_radius=BOX_CORNER_RADIUS):
    """Draw a rectangle with rounded corners"""
    if pygame.version.vernum[0] >= 2:
        # Use built-in rounded rect in Pygame 2.0+
        pygame.draw.rect(screen, color, rect, border_radius=border_radius)
    else:
        # Fallback for older Pygame versions
        pygame.draw.rect(screen, color, rect)
