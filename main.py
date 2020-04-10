from typing import List, Tuple, TextIO

import os.path
import sys

from scipy.spatial import Delaunay
import numpy as np

from classes import Point, VoronoyPoint, VoronoyEdge
from graph_utils import Graph, create_graph_from_edges, compute_tunnel
from math import radians, cos, sin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d

TunnelCenterLine: type = List[VoronoyPoint]


def load_data(filename: str) -> List[Point]:
    """
    Чтение файла
    :param filename: путь к файлу типа pdb
    :return: список координат атомов в этом файле
    """
    file: TextIO = open(filename)
    res = []
    while file.readable():
        line: str = file.readline()
        if line.startswith("ENDMDL") or line == "":
            break
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if line[-3] == 'H' or line[-4] == 'H':  # водород
            continue

        point_coords = line[31:54].split()
        x = float(point_coords[0])
        y = float(point_coords[1])
        z = float(point_coords[2])

        point = Point(x, y, z)
        res.append(point)

    return res


def get_centroid(tetrahedron) -> VoronoyPoint:
    """
    Поиск центроида тетраэдра
    :param tetrahedron: тетраэдр, заданный как массив точек (массивов координат)
    :return: Точка с координатами и списком соседних точек
    """
    tetrahedron_vertices: List[Point] = list(map(lambda arr: Point(arr[0], arr[1], arr[2]), tetrahedron))

    x: float = sum(map(lambda point: point.coords[0], tetrahedron_vertices)) / 4
    y: float = sum(map(lambda point: point.coords[1], tetrahedron_vertices)) / 4
    z: float = sum(map(lambda point: point.coords[2], tetrahedron_vertices)) / 4
    return VoronoyPoint(x, y, z, tetrahedron_vertices)


def get_voronoi_vertex(tetrahedron: List[List[float]]) -> VoronoyPoint:
    tetrahedron_vertices: List[Point] = list(map(lambda arr: Point(arr[0], arr[1], arr[2]), tetrahedron))

    a = tetrahedron[0]
    b = tetrahedron[1]
    c = tetrahedron[2]
    s = tetrahedron[3]

    xs = s[0]
    ys = s[1]
    zs = s[2]

    xa = a[0]
    ya = a[1]
    za = a[2]

    xb = b[0]
    yb = b[1]
    zb = b[2]

    xc = c[0]
    yc = c[1]
    zc = c[2]

    sa = xs * xs - xa * xa + ys * ys - ya * ya + zs * zs - za * za
    sb = xs * xs - xb * xb + ys * ys - yb * yb + zs * zs - zb * zb
    sc = xs * xs - xc * xc + ys * ys - yc * yc + zs * zs - zc * zc

    x_mat = [[sa, ys - ya, zs - za],
             [sb, ys - yb, zs - zb],
             [sc, ys - yc, zs - zc],
             ]

    y_mat = [[xs - xa, sa, zs - za],
             [xs - xb, sb, zs - zb],
             [xs - xc, sc, zs - zc],
             ]

    z_mat = [[xs - xa, ys - ya, sa],
             [xs - xb, ys - yb, sb],
             [xs - xc, ys - yc, sc],
             ]

    del_mat = [
        [xs - xa, ys - ya, zs - za],
        [xs - xb, ys - yb, zs - zb],
        [xs - xc, ys - yc, zs - zc],
    ]

    znam = 2 * np.linalg.det(del_mat)

    if abs(znam) > 1e-6 and abs(znam) < 1e4:
        x0 = np.linalg.det(x_mat) / znam
        y0 = np.linalg.det(y_mat) / znam
        z0 = np.linalg.det(z_mat) / znam
    else:
        x0 = (xs + xa + xb + xc) / 4
        y0 = (ys + ya + yb + yc) / 4
        z0 = (zs + za + zb + zc) / 4

    return VoronoyPoint(x0, y0, z0, tetrahedron_vertices)


def plot(points: List[Point], tri: Delaunay, vertices: List[VoronoyPoint], edges: List[VoronoyEdge]):
    fig = plt.figure()
    # fig.set_tight_layout(True)
    ax = Axes3D(fig)  # fig.add_subplot(111, projection='3d')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    x_points = [p.coords[0] for p in points]
    y_points = [p.coords[1] for p in points]
    z_points = [p.coords[2] for p in points]

    # визуализация триангляции
    # ax.plot_trisurf(x_points, y_points, z_points, triangles=tri.simplices, cmap=plt.cm.Spectral)
    for vertex in vertices:
        neibs = vertex.neighbours

        tmp = [*neibs, neibs[0], neibs[2], neibs[1], neibs[3]]

        x_points_ = [p.coords[0] for p in tmp]
        y_points_ = [p.coords[1] for p in tmp]
        z_points_ = [p.coords[2] for p in tmp]

        tri = art3d.Poly3DCollection([neib.coords for neib in neibs])
        tri.set_alpha(0.1)
        # tri.set_color('grey')
        # ax.add_collection3d(tri)

        ax.plot(x_points_, y_points_, z_points_, alpha=0.1)

    edge_mids = [(edge.start_point + edge.finish_point) / 2 for edge in edges]
    x_points_e = [p.coords[0] for p in edge_mids]
    y_points_e = [p.coords[1] for p in edge_mids]
    z_points_e = [p.coords[2] for p in edge_mids]
    ax.scatter(x_points_e, y_points_e, z_points_e, marker='*', color='red')
    # исходные вершины
    ax.scatter(x_points, y_points, z_points, marker='*', color='green')

    # вершины Вороного(не shallow)
    x_points = [vertices[i].coords[0] for i in range(len(vertices))]
    y_points = [vertices[i].coords[1] for i in range(len(vertices))]
    z_points = [vertices[i].coords[2] for i in range(len(vertices))]

    ax.scatter(x_points, y_points, z_points, marker='v', color='blue')

    # for tunnel in tunnels:
    #     x = [point.coords[0] for point in tunnel]
    #     y = [point.coords[1] for point in tunnel]
    #     z = [point.coords[2] for point in tunnel]
    #     ax.plot(x, y, z, color='yellow')

    plt.show()


def get_voronoy_edges(tri: Delaunay) -> Tuple[List[VoronoyPoint], List[VoronoyEdge]]:
    """
    Ребра для Графа Вороного
    :param tri: результаты триангулция Делануа
    :return: Граф Вороного - кортеж из списка вершин и списка ребер
    """
    p = tri.points[tri.vertices]
    centroids: List[VoronoyPoint] = list(map(get_voronoi_vertex, p))
    lines: List[VoronoyEdge] = []
    for index, neighbours_index in enumerate(tri.neighbors):
        for neighbour in neighbours_index:
            if neighbour > index:  # без повторов
                lines.append(VoronoyEdge(centroids[index], centroids[neighbour], index, neighbour))
    return centroids, lines


def dump2pdb(pdb_file: TextIO, data: TunnelCenterLine) -> None:
    """
    Сохранение в pdb линии туннеля
    :param pdb_file: открытый поток файла
    :param data: туннел заданный списком точек
    :return: void
    """
    pdb_fmt = "{:6s}{:5d}  {:<4s}{:>3s} {:1s}{:4d}    {:>8s}{:>8s}{:>8s}                {:>8s}"
    fmtc = "{:8.3f}"
    for ia, r in enumerate(data):
        c1 = fmtc.format(r.coords[0])
        c2 = fmtc.format(r.coords[1])
        c3 = fmtc.format(r.coords[2])
        ch = " "
        pdb_file.write(pdb_fmt.format("HETATM", ia % 100000, 'X', "DUM", ch, ia % 10000, c1, c2, c3, 'I'))
        pdb_file.write("\n")


def find_center_line(input_file: str, s_start: Point, min_z: float, max_z: float) -> TunnelCenterLine:
    """
    Поиск туннеля в молекуле
    :param input_file: название файла с туннелем(должен быть в папке input_files)
    :param s_start: начальная точка x, y, z
    :param min_z: нижняя граница туннеля
    :param max_z: верхняя граница туннеля
    :return:туннель заданный последовательностью точек
    """
    directory: str = "input_files"
    points: List[Point] = load_data(os.path.join(directory, input_file))

    # r: float = 100
    # circle_ponts = [(r * cos(radians(teta + 45)), r * sin(radians(teta + 45))) for teta in range(360) if teta % 2 == 0]
    # points: List[Point] = [Point(x, y, z) for x, y in circle_ponts for z in range(-50, 51, 10)]

    save_center_line_to_pdb(points, 'in.pdb')

    np_points = np.array(list(map(lambda p: p.coords, points)))

    tri = Delaunay(np_points)

    vertices, edges = get_voronoy_edges(tri)
    g: int = len(vertices)
    save_center_line_to_pdb(vertices, 'verts.pdb')

    graph: Graph = create_graph_from_edges(edges, set(range(g)), set(range(len(edges))))

    s_start_ind: int = min(range(g), key=lambda i: s_start.dist(vertices[i]))

    # def stop(ind: int) -> bool:
    #     return vertices[ind].coords[2] > max_z or vertices[ind].coords[2] < min_z

    stop_max = lambda ind: vertices[ind].coords[2] > max_z
    stop_min = lambda ind: vertices[ind].coords[2] < min_z
    best_tunnel_indexes: List[int] = compute_tunnel(graph, g, s_start_ind, stop_min)
    #best_tunnel: List[Point] = [vertices[i] for i in best_tunnel_indexes]

    best_tunnel_indexes2: List[int] = compute_tunnel(graph, g, s_start_ind, stop_max)
    #best_tunnel2: List[Point] = [vertices[i] for i in best_tunnel_indexes2]
    tunnel_indexes = [*best_tunnel_indexes, *(list(reversed(best_tunnel_indexes2))[1:])]
    best_tunnel:List[Point] = [vertices[tunnel_indexes[0]]]
    for i in range(1, len(tunnel_indexes)):
        mid_point = (vertices[tunnel_indexes[i-1]]+vertices[tunnel_indexes[i]])/2
        best_tunnel.append(mid_point)
        best_tunnel.append(vertices[tunnel_indexes[i]])
    #tunnel = [*best_tunnel, *best_tunnel2]
    #plot(points, tri, vertices, edges)
    return best_tunnel


def save_center_line_to_pdb(center_line: TunnelCenterLine, filename: str = "out.pdb") -> None:
    """
    Сохранение туннеля в файл
    :param center_line: туннель
    :param filename: имя файла
    :return: void
    """
    file: TextIO = open(filename, 'w')
    file.write("MODEL 1\n")
    dump2pdb(file, center_line)
    file.write("ENDMODEL")
    file.close()


def main():
    if len(sys.argv) != 7:
        print('Ошибка, необходимо 6 параметров')
        return
    # TODO обработка ошибок
    filename: str = sys.argv[1]

    start_x: float = float(sys.argv[2])
    start_y: float = float(sys.argv[3])
    start_z: float = float(sys.argv[4])

    min_z: float = float(sys.argv[5])
    max_z: float = float(sys.argv[6])
    s_start = Point(start_x, start_y, start_z)

    center_line: TunnelCenterLine = find_center_line(filename, s_start, min_z, max_z)
    save_center_line_to_pdb(center_line)


if __name__ == '__main__':
    # TODO убрать профилировщик
    import cProfile

    cProfile.run("main()")
