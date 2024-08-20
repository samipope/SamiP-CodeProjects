from __future__ import annotations
import math as math
__all__ = ['CanvasBase', 'math']
class CanvasBase:
    """
    Base class for specialized canvas backends
    """
    def _getLinePoints(self, p1, p2, dash):
        ...
    def addCanvasDashedWedge(self, p1, p2, p3, dash = (2, 2), color = (0, 0, 0), color2 = None, **kwargs):
        """
        Draw a dashed wedge
        
                   The wedge is identified by the three points `p1`, `p2`, and `p3`.
                   It will be drawn using the given `color`; if `color2` is specified
                   it will be used for the second half of the wedge
        
                   TODO: fix comment, I'm not sure what `dash` does
        
                
        """
    def addCanvasLine(self, p1, p2, color = (0, 0, 0), color2 = None, **kwargs):
        """
        Draw a single line on the canvas
        
                    This function will draw a line between `p1` and `p2` with the
                    given `color`.
                    If `color2` is specified, it will be used to draw the second half
                    of the segment
                
        """
    def addCanvasPolygon(self, ps, color = (0, 0, 0), **kwargs):
        """
        Draw a polygon
        
                   Draw a polygon identified by vertexes given in `ps` using
                   the given `color`
                
        """
    def addCanvasText(self, text, pos, font, color = (0, 0, 0), **kwargs):
        """
        Draw some text
        
                   The provided `text` is drawn at position `pos` using the given
                   `font` and the chosen `color`.
                
        """
    def flush(self):
        """
        Complete any remaining draw operation
        
                   This is supposed to be the last operation on the canvas before
                   saving it
                
        """
