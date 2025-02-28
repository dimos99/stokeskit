import pygame
import time
from config import (BUTTON_BG_COLOR, BUTTON_HOVER_COLOR, BUTTON_ACTIVE_COLOR,
                   BUTTON_TEXT_COLOR, BUTTON_BORDER_COLOR, BUTTON_CORNER_RADIUS,
                   TEXT_COLOR, TEXT_BG_COLOR)

class Button:
    def __init__(self, x, y, width, height, text, toggle=False, active=False, tooltip=None):
        self.rect = pygame.Rect(x, y, width, height)
        self.text = text
        self.is_toggle = toggle
        self.active = active if toggle else False
        self.hover = False
        self.tooltip = tooltip
        self.font = pygame.font.SysFont('Arial', 14, True)
        self.tooltip_font = pygame.font.SysFont('Arial', 12)
        self.tooltip_shown = False
        self.tooltip_timer = 0
        
    def draw(self, surface):
        # Determine button color based on state
        if self.is_toggle and self.active:
            color = BUTTON_ACTIVE_COLOR
        elif self.hover:
            color = BUTTON_HOVER_COLOR
        else:
            color = BUTTON_BG_COLOR
            
        # Draw button background
        pygame.draw.rect(surface, color, self.rect, border_radius=BUTTON_CORNER_RADIUS)
        pygame.draw.rect(surface, BUTTON_BORDER_COLOR, self.rect, 1, border_radius=BUTTON_CORNER_RADIUS)
        
        # Draw button text
        text_surface = self.font.render(self.text, True, BUTTON_TEXT_COLOR)
        text_rect = text_surface.get_rect(center=self.rect.center)
        surface.blit(text_surface, text_rect)
        
        # Draw tooltip if hovering
        if self.hover and self.tooltip and self.tooltip_shown:
            tooltip_surface = self.tooltip_font.render(self.tooltip, True, TEXT_COLOR)
            tooltip_rect = tooltip_surface.get_rect()
            tooltip_rect.topleft = (self.rect.right + 5, self.rect.top)
            
            # Draw tooltip background
            bg_rect = tooltip_rect.inflate(10, 6)
            pygame.draw.rect(surface, TEXT_BG_COLOR, bg_rect, border_radius=3)
            
            # Draw tooltip text
            surface.blit(tooltip_surface, tooltip_rect)
        
    def update(self, mouse_pos, mouse_clicked=False):
        prev_hover = self.hover
        self.hover = self.rect.collidepoint(mouse_pos)
        
        # Start tooltip timer when starting hover
        if self.hover and not prev_hover:
            self.tooltip_timer = time.time()
            self.tooltip_shown = False
        
        # Show tooltip after hovering for 0.7 seconds
        if self.hover and not self.tooltip_shown and (time.time() - self.tooltip_timer > 0.7):
            self.tooltip_shown = True
        
        # Hide tooltip when not hovering
        if not self.hover:
            self.tooltip_shown = False
        
        # Handle click
        if mouse_clicked and self.hover:
            if self.is_toggle:
                self.active = not self.active
            return True
        return False
