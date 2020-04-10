from typing import Tuple
from math import sqrt


class Point:
    """
    3D точка
    coords = [x, y, z]
    """
    coords: Tuple[float, float, float]

    def __repr__(self):
        return f"Point {self.coords[0]}, {self.coords[1]}, {self.coords[2]}"

    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        self.coords = x, y, z

    def dist(self, other: "Point") -> float:
        a = self.coords
        b = other.coords
        tmp = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
        return sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2])

    def __truediv__(self, other: float) -> "Point":
        coords = [self.coords[0] / other, self.coords[1] / other, self.coords[2] / other]
        return Point(coords[0], coords[1], coords[2])

    def __rmul__(self, other: float):
        coords = [self.coords[0] * other, self.coords[1] * other, self.coords[2] * other]
        return Point(coords[0], coords[1], coords[2])

    def __sub__(self, other: "Point") -> "Point":
        tmp = [self.coords[0] - other.coords[0], self.coords[1] - other.coords[1], self.coords[2] - other.coords[2]]
        return Point(tmp[0], tmp[1], tmp[2])

    def __add__(self, other: "Point") -> "Point":
        coords = [self.coords[0] + other.coords[0], self.coords[1] + other.coords[1], self.coords[2] + other.coords[2]]
        return Point(coords[0], coords[1], coords[2])

    def __abs__(self) -> float:
        tmp = [self.coords[0] + self.coords[0], self.coords[1] + self.coords[1], self.coords[2] + self.coords[2]]
        return sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2])

    def __hash__(self):
        return hash(self.coords)

    def __eq__(self, other):
        return all(self.coords[i] == other.coords[i] for i in range(3))

    def copy(self):
        return Point(self.coords[0], self.coords[1], self.coords[2])
