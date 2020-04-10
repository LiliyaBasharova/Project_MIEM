from typing import List
from .point import Point


class VoronoyPoint(Point):
    """
    Вершина в графе Вороного
    coords = [x, y, z]

    """
    neighbours: List[Point]

    def __init__(self, x: float, y: float, z: float, neibhours: List[Point]):
        super(VoronoyPoint, self).__init__(x, y, z)
        self.neighbours = neibhours

    def __repr__(self):
        return f"VoronoiPoint {self.coords[0]}, {self.coords[1]}, {self.coords[2]}"
