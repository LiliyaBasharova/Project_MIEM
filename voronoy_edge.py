from typing import List
from .point import Point
from .voronoy_vertex import VoronoyPoint


class Edge(object):
    """
    Ребро графа
    start_index - индекс начальной вершины
    finish_index - индекс конечной вершины
    cost - стоимость
    """
    start_index: int
    finish_index: int
    cost: float

    def __init__(self, start: int, finish: int):
        self.start_index = start
        self.finish_index = finish

    def __repr__(self):
        return f"Edge from {self.start_index} to {self.finish_index} cost = {self.cost}"


class VoronoyEdge(Edge):
    """
    Ребро в графе Вороного
    start_point - начальная вершина
    finish_point - конечная вершина
    """
    start_point: VoronoyPoint
    finish_point: VoronoyPoint
    max_radius_sphere: float

    def __init__(self, start: VoronoyPoint, finish: VoronoyPoint, st_index: int, fi_index: int):
        super().__init__(st_index, fi_index)
        self.start_point = start
        self.finish_point = finish
        self.cost = self.edge_count()

    def edge_count(self) -> float:
        """
        Расчёт стоимости ребра(узкое место!!)
        :return:
        """
        start_copy: Point = self.start_point.copy()
        finish_copy: Point = self.finish_point.copy()

        all_points: List[Point] = list(set(self.start_point.neighbours) & set(self.finish_point.neighbours))
        tetrahedron_vertex = all_points[0]
        vec: Point = finish_copy - start_copy

        z: float = 1.5
        N: int = 100
        step = abs(vec) / N

        integral: float = 0
        self.max_radius_sphere = 1000000
        for _ in range(N):
            start_copy += step * vec
            max_r = tetrahedron_vertex.dist(start_copy) #min([start_copy.dist(atom) for atom in all_points])
            self.max_radius_sphere = min(self.max_radius_sphere, max_r)
            integral += pow(max_r, -z) * step

        return integral
