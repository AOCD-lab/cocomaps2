"""
Defines the Ring_params class, which computes and stores geometric parameters for
aromatic ring systems detected in molecular structures. The class calculates each
ring’s centroid, plane coefficients, and polygonal boundaries based on atom coordinates,
allowing accurate geometric reasoning about π-related interactions.

It also provides methods to project nearby atoms or points onto ring planes and
evaluate whether those points lie within or near the π-system boundary, based on
distance thresholds specific to interaction types such as cation–π, anion–π, or
lone pair–π contacts. This enables the geometric validation of π-centered
noncovalent interactions in complex biomolecular assemblies.
"""

from typing import Optional
from utils import (
    create_plane,
    project_point_onto_plane,
    get_centroid,
    is_point_inside_3d_polygon,
)
from update_constants import tc
from constants import (
    CHAIN1_IDENTIFIER,
    CHAIN2_IDENTIFIER,
    LONE_PAIR_PI_PRINT_NAME,
    ANION_PI_PRINT_NAME,
    CATION_PI_PRINT_NAME,
    POLAR_PRINT_NAME,
)
from utils import get_dist_between_points


ANION_PI_DIST = tc.ANION_PI_DIST
CATION_PI_DIST = tc.CATION_PI_DIST
LONEPAIR_PI_DIST = tc.LONEPAIR_PI_DIST


class Ring_params:
    def __init__(self) -> None:
        pass

    def make_chains_ring_comps(self, rings: dict):
        self.polygons_coords = {}
        for chain, res_num_dicts in rings.items():
            if chain not in self.polygons_coords.keys():
                self.polygons_coords[chain] = {}
            for res_num, atoms in res_num_dicts.items():
                cuur_ring_coords = []
                for atom, coords in atoms.items():
                    cuur_ring_coords.append(coords)
                self.polygons_coords[chain][res_num] = cuur_ring_coords
                if [] in cuur_ring_coords:
                    del self.polygons_coords[chain][res_num]

        self.centroids = {}
        self.planes_coefficients = {}
        for chain, res_num_dicts in self.polygons_coords.items():
            if chain not in self.centroids:
                self.centroids[chain] = {}
            if chain not in self.planes_coefficients:
                self.planes_coefficients[chain] = {}
            for res_num, polygon_coords in res_num_dicts.items():
                temp_centroid = get_centroid(polygon_coords)
                self.centroids[chain][res_num] = temp_centroid
                self.planes_coefficients[chain][res_num] = create_plane(
                    temp_centroid, polygon_coords[0], polygon_coords[1]
                )
                temp_centroid = None

    def check_point_inside(
        self, point_to_project, plane_coefficients, polygon_coords, type_, polygon_centroid:Optional[list] = None
    ):
        distance_return = None
        projection, distance1 = project_point_onto_plane(
            point_to_project, plane_coefficients
        )
        threshold = 0
        # distance = get_dist_between_points(point_to_project, polygon_centroid)
        distance1 = round(distance1,2)
        if type_ == LONE_PAIR_PI_PRINT_NAME:
            threshold = LONEPAIR_PI_DIST
        elif type_ == ANION_PI_PRINT_NAME:
            threshold = ANION_PI_DIST
        elif type_ == CATION_PI_PRINT_NAME:
            threshold = CATION_PI_DIST
        else:
            threshold = 4
        if distance1 <= threshold:
            inside_bool = is_point_inside_3d_polygon(projection, polygon_coords)
            if inside_bool:
                distance_return = distance1
            return distance_return
        else:
            return False
