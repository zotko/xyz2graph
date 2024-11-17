"""Configuration and factory classes for Plotly visualization buttons.

This module provides classes and functions for creating and managing buttons
used in molecular structure visualization controls.

Classes:
    ButtonConfig: Configuration class for button positioning and styling.
    ButtonFactory: Factory class for creating different types of Plotly button configurations.
"""

from typing import Any, Dict, List, TypedDict

from .layout_config import LAYOUT


class ButtonConfig(TypedDict, total=False):
    """Configuration dictionary for button appearance and behavior.

    Attributes:
        direction: Direction of button layout ("right", "left", "up", "down")
        pad_r: Right padding in pixels
        pad_t: Top padding in pixels
        showactive: Whether to highlight the active button
        x: X-position of button group (0-1)
        y: Y-position of button group (0-1)
        xanchor: X-anchor point ("left", "center", "right")
        yanchor: Y-anchor point ("top", "middle", "bottom")
    """

    direction: str
    pad_r: int
    pad_t: int
    showactive: bool
    x: float
    y: float
    xanchor: str
    yanchor: str


DEFAULT_BUTTON_CONFIG = ButtonConfig(
    direction="right",
    pad_r=10,
    pad_t=10,
    showactive=True,
    x=LAYOUT["common"]["button_x"],
    y=LAYOUT["background"]["button_y"],  # Use layout settings for default position
    xanchor="left",
    yanchor="top",
)


class ButtonFactory:
    """Factory class for creating Plotly button configurations."""

    @staticmethod
    def create_toggle_button_group(
        title: str,
        y_position: float,
        buttons_config: List[Dict[str, Any]],
        active_index: int = 0,
        button_config: ButtonConfig = DEFAULT_BUTTON_CONFIG,
    ) -> tuple[Dict[str, Any], Dict[str, Any]]:
        """Create a button group with title annotation.

        Args:
            title: Title text for the button group
            y_position: Y-position for the button group (0-1)
            buttons_config: List of button configurations
            active_index: Index of initially active button
            button_config: Configuration for button appearance and behavior

        Returns:
            Tuple containing:
                - Button group configuration dictionary
                - Title annotation configuration dictionary
        """
        button_group = {
            "buttons": buttons_config,
            "direction": button_config["direction"],
            "pad": {"r": button_config["pad_r"], "t": button_config["pad_t"]},
            "showactive": button_config["showactive"],
            "active": active_index,
            "x": button_config["x"],
            "y": y_position,
            "xanchor": button_config["xanchor"],
            "yanchor": button_config["yanchor"],
            "type": "buttons",
        }

        title_annotation = {
            "text": title,
            "x": LAYOUT["common"]["text_x"],
            "xref": "paper",
            "y": y_position - 0.02,  # Offset title slightly below buttons
            "yref": "paper",
            "align": "left",
            "showarrow": False,
        }

        return button_group, title_annotation

    @staticmethod
    def create_background_toggle() -> tuple[Dict[str, Any], Dict[str, Any]]:
        """Create a background visibility toggle button group.

        Returns:
            Tuple containing:
                - Background toggle button configuration
                - Title annotation configuration
        """
        buttons_config = [
            {
                "label": "Show",
                "method": "relayout",
                "args": [
                    {
                        "scene.xaxis.showbackground": True,
                        "scene.yaxis.showbackground": True,
                        "scene.zaxis.showbackground": True,
                    }
                ],
            },
            {
                "label": "Hide",
                "method": "relayout",
                "args": [
                    {
                        "scene.xaxis.showbackground": False,
                        "scene.yaxis.showbackground": False,
                        "scene.zaxis.showbackground": False,
                    }
                ],
            },
        ]

        return ButtonFactory.create_toggle_button_group(
            title=LAYOUT["background"]["title"],
            y_position=LAYOUT["background"]["button_y"],
            buttons_config=buttons_config,
        )
