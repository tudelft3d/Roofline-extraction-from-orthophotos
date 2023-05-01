import numpy as np
import copy
import os
os.environ.setdefault('OPENCV_IO_MAX_IMAGE_PIXELS', pow(2,50).__str__())
import cv2

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from shapely.geometry import Polygon as ShapeL_polygon
from shapely.geometry import MultiPolygon as ShapeL_multipolygon
from shapely.geometry import Point as ShapeL_point
from shapely.ops import unary_union as ShapeL_unary_union
from shapely.ops import polygonize as ShapeL_polygonize
import rtree

class WPolygon():
    def __init__(self, id, building_id):
        self.id = id
        self.building_id = building_id
        self.bbox_coord = None
        self.geo_points = []
        self.pix_points = []
        self.geo_center = None
        self.pix_center = None
        self.buffered_pix_points = []
        self.edges = []

    def convert_geo_pt_to_pix_pt(self, pts, top_left, pix_size, bbox_pix):
        pts_pix = (int((pts[0] - top_left[0]) / pix_size) - bbox_pix[0], int((top_left[1] - pts[1]) / pix_size) - bbox_pix[1])
        return pts_pix

    def convert_from_geo_to_pix(self, top_left, pix_size, bbox_pix):
        self.pix_points = []
        if len(self.geo_points) != 0:
            for pts in self.geo_points:
                pts_pix = (int((pts[0] - top_left[0]) / pix_size) - bbox_pix[0], int((top_left[1] - pts[1]) / pix_size) - bbox_pix[1])
                self.pix_points.append(pts_pix)
            self.pix_points = np.array(self.pix_points)
        else:
            self.pix_points = np.zeros((0, 2), dtype=np.float64)

    def convert_from_pix_to_geo(self, top_left, pix_size, bbox_pix):
        self.geo_points = []
        if len(self.pix_points) != 0:
            for pts_pix in self.pix_points:
                pts = (top_left[0] + float((pts_pix[0] + bbox_pix[0]) * pix_size), top_left[1] - float((pts_pix[1] + bbox_pix[1]) * pix_size))
                self.geo_points.append(pts)
            self.geo_points = np.array(self.geo_points)
        else:
            self.geo_points = np.zeros((0, 2), dtype=np.float64)

    def buffering_polygon(self, args):
        shp_poly = ShapeL_polygon(self.pix_points)
        # extend the polygon by a certain value
        extshp_poly = shp_poly.buffer(args.polygon_buffer)
        self.buffered_pix_points = [(int(x), int(y)) for x, y in extshp_poly.exterior.coords]

    def draw_polygon(self, img, draw_buffer):
        used_points = self.pix_points
        if draw_buffer > 1:
            used_points = self.buffered_pix_points
        if len(used_points) != 0:
            pi = 0
            for pts_pix in used_points:
                if pi < len(used_points) - 1:
                    img = cv2.circle(img, pts_pix, 2, (0, 255, 0), 0)
                    img = cv2.line(img, pts_pix, used_points[pi + 1], (0, 0, 255), 1)
                pi += 1

class WPolygons():
    def __init__(self):
        self.file = None
        self.shapes = None #shapes read from input
        self.layer_line = None
        self.orig_building_polygon_dict = dict()
        self.orig_building_shape_dict = dict()
        self.building_polygon_dict = dict()
        self.building_bbox_dict = dict()
        self.building_gt_rooflines_dict = dict()
        self.merged_to_original_map = None
        self.srs = None

    def read_poly(self, path, args):
        self.read_gpkg(path)
        self.preprocess_raw_gpkg(args)

    def initialize_gpkg_for_write(self, path, building_polys):
        path = path + ".gpkg"
        driver = ogr.GetDriverByName("GPKG")
        self.file = driver.CreateDataSource(path)
        self.layer_line = self.file.CreateLayer("lines", building_polys.srs, ogr.wkbLineString)
        # Add file name to the new layer
        building_polys_layer = building_polys.file.GetLayer()
        building_polys_layer_def = building_polys_layer.GetLayerDefn()
        for i in range(building_polys_layer_def.GetFieldCount()):
            building_polys_field_def = building_polys_layer_def.GetFieldDefn(i)
            building_polys_field_name = building_polys_field_def.GetName()
            if building_polys_field_name == 'fid':
                building_polys_field_name = 'fid_pand'
            building_polys_field_type = building_polys_field_def.GetType()
            field_new = ogr.FieldDefn(building_polys_field_name, building_polys_field_type)
            self.layer_line.CreateField(field_new)

    def read_gpkg(self, path):
        self.file = ogr.Open(path, gdal.GA_ReadOnly)
        # Check if the file was successfully opened
        if self.file is None:
            print("Error opening the GeoPackage file:", gdal.GetLastErrorMsg())
        # else:
        #     print("GeoPackage file successfully opened")

    def close_file(self, args):
        self.file = None
        self.layer_line = None

    def add_lines(self, id, poly, srs):
        layer_line = self.file.CreateLayer("lines", srs, ogr.wkbLineString)
        # Add a new line feature to the layer
        line = ogr.Geometry(ogr.wkbLineString)
        for ed in poly.edges:
            line.AddPoint(poly.geo_points[ed[0]], poly.geo_points[ed[1]])
        feature = ogr.Feature(layer_line.GetLayerDefn())
        feature.SetGeometry(line)
        layer_line.CreateFeature(feature)

    def add_partition_lines(self, id, sing_b, used_shapes):
        # Add a new line feature to the layer
        ei = 0
        for ed in sing_b.edges:
            pt1 = sing_b.geo_points[ed[0]]
            pt2 = sing_b.geo_points[ed[1]]

            # Create a new feature
            featureDefn = self.layer_line.GetLayerDefn()
            feature = ogr.Feature(featureDefn)

            # Set the geometry of the feature to a new line
            line = ogr.Geometry(ogr.wkbLineString)
            line.AddPoint(pt1[0], pt1[1])
            line.AddPoint(pt2[0], pt2[1])
            feature.SetGeometry(line)

            # Add the feature to the layer
            self.layer_line.CreateFeature(feature)

            # Add attributes
            old_feature_ei = used_shapes[sing_b.edge_shape_ind[ei]]
            for i in range(old_feature_ei.GetFieldCount()):
                #Get old file details
                old_field_def_ei = old_feature_ei.GetFieldDefnRef(i)
                old_field_name_ei = old_field_def_ei.GetName()
                old_field_value_ei = old_feature_ei.GetField(old_field_name_ei)
                if old_field_name_ei == 'fid':
                    old_field_name_ei = 'fid_pand'
                feature.SetField(old_field_name_ei, old_field_value_ei)
            self.layer_line.SetFeature(feature)
            # Cleanup
            feature.Destroy()
            ei += 1
    def get_raw_polygons_for_merge(self, building_id, layer_name):
        id = 0
        polys_b_id = []
        polys_geo_points = []
        for shape in self.shapes:
            b_id = shape.GetField(building_id)
            polys_b_id.append(b_id)
            poly = WPolygon(id, b_id)
            building_geo = shape.GetGeometryRef()
            if building_geo.GetGeometryType() == ogr.wkbPolygon or \
                    building_geo.GetGeometryName() == "MULTIPOLYGON":
                ring = building_geo.GetGeometryRef(0)
                num_points = ring.GetPointCount()
                poly.geo_points = np.empty((num_points, 2), dtype=np.float64)
                x_avg = 0
                y_avg = 0
                for i in range(num_points):
                    x, y, z = ring.GetPoint(i)
                    poly.geo_points[i, 0] = x
                    poly.geo_points[i, 1] = y
                    x_avg += x
                    y_avg += y
                polys_geo_points.append(poly.geo_points)
                x_avg /= num_points
                y_avg /= num_points
                poly.geo_center = (x_avg, y_avg)
                if b_id not in self.orig_building_polygon_dict.keys():
                    self.orig_building_polygon_dict[b_id] = [poly]
                    self.orig_building_shape_dict[b_id] = [shape]
                else:
                    self.orig_building_polygon_dict[b_id].append(poly)
                    self.orig_building_shape_dict[b_id].append(shape)
                id += 1
        return polys_b_id, polys_geo_points

    def get_polygon_details(self, poly):
        # Get the bounding box of the polygon
        bbox = poly.bounds
        # Get the ordered points of the polygon
        points = np.array(poly.exterior.coords)
        return bbox, points
    def merge_adjacent_polygons(self, polys_b_id, polys_geo_points, building_buffer_dis):
        # Convert each numpy array to a shapely Polygon object
        polygons = [ShapeL_polygon(poly) for poly in polys_geo_points]
        # Create a MultiPolygon object
        mpoly = ShapeL_multipolygon(polygons)
        # Use buffer() to ensure adjacent polygons share the same geometry
        buffer_distance = building_buffer_dis
        merged_poly = mpoly.buffer(buffer_distance)
        # Merge all polygons into a single polygon
        merged_poly = ShapeL_unary_union(merged_poly)
        # Polygonize the merged polygon to get a list of polygons
        merged_polys = list(ShapeL_polygonize(merged_poly))
        # Iterate over the merged polygons and assign the first id, and create a map of ids between the original polygons
        details = []
        new_ids = []
        merged_to_original = {}
        poly_index = 0
        for poly_m in merged_polys:
            valid_poly = False
            first_poly_add = False
            for j, poly_o in enumerate(polygons):
                if poly_o.centroid.within(poly_m):
                    valid_poly = True
                    if poly_index not in merged_to_original:
                        merged_to_original[poly_index] = []
                    merged_to_original[poly_index].append(polys_b_id[j])
                    if not first_poly_add:
                        first_poly_add = True
                        new_ids.append(polys_b_id[j])
            if valid_poly:
                details.append(self.get_polygon_details(poly_m))
                poly_index += 1

        # Extract the bounding box and points of each polygon
        bboxes, points_list = zip(*details)
        return new_ids, bboxes, points_list, merged_to_original

    def process_to_merge_adjacent_polygons(self, building_id, layer_name, building_buffer_dis):
        polys_b_id, polys_geo_points = self.get_raw_polygons_for_merge(building_id, layer_name)
        polys_b_id, bboxes, poly_points_list, self.merged_to_original_map = self.merge_adjacent_polygons(polys_b_id, polys_geo_points, building_buffer_dis)
        id = 0
        # Assign ids to merged polygons
        for b_id, box, pts in zip(polys_b_id, bboxes, poly_points_list):
            poly = WPolygon(id, b_id)
            poly.bbox_coord = list(box)  # (lower left, upper right)
            poly.geo_points = pts
            if b_id not in self.building_polygon_dict.keys():
                self.building_polygon_dict[b_id] = [poly]
                self.building_bbox_dict[b_id] = poly.bbox_coord
            else:
                self.building_polygon_dict[b_id].append(poly)
                if poly.bbox_coord[0] < self.building_bbox_dict[b_id][0]:
                    self.building_bbox_dict[b_id][0] = poly.bbox_coord[0]
                if poly.bbox_coord[1] < self.building_bbox_dict[b_id][1]:
                    self.building_bbox_dict[b_id][1] = poly.bbox_coord[1]
                if poly.bbox_coord[2] > self.building_bbox_dict[b_id][2]:
                    self.building_bbox_dict[b_id][2] = poly.bbox_coord[2]
                if poly.bbox_coord[3] > self.building_bbox_dict[b_id][3]:
                    self.building_bbox_dict[b_id][3] = poly.bbox_coord[3]
            id += 1

    def process_raw_shapes(self, building_id, layer_name):
        id = 0
        for shape in self.shapes:
            b_id = shape.GetField(building_id)
            poly = WPolygon(id, b_id)
            building_geo = shape.GetGeometryRef()
            if building_geo.GetGeometryType() == ogr.wkbPolygon or \
                    building_geo.GetGeometryName() == "MULTIPOLYGON":
                xmin, xmax, ymin, ymax = building_geo.GetEnvelope()
                poly.bbox_coord = [xmin, ymin, xmax, ymax]  # (lower left, upper right)
                ring = building_geo.GetGeometryRef(0)
                num_points = ring.GetPointCount()
                poly.geo_points = np.empty((num_points, 2), dtype=np.float64)
                x_avg = 0
                y_avg = 0
                for i in range(num_points):
                    x, y, z = ring.GetPoint(i)
                    poly.geo_points[i, 0] = x
                    poly.geo_points[i, 1] = y
                    x_avg += x
                    y_avg += y
                x_avg /= num_points
                y_avg /= num_points
                poly.geo_center = (x_avg, y_avg)
            if b_id not in self.building_polygon_dict.keys():
                self.building_polygon_dict[b_id] = [poly]
                self.orig_building_shape_dict[b_id] = [shape]
                self.building_bbox_dict[b_id] = poly.bbox_coord
            else:
                self.building_polygon_dict[b_id].append(poly)
                self.orig_building_shape_dict[b_id].append(shape)
                if poly.bbox_coord[0] < self.building_bbox_dict[b_id][0]:
                    self.building_bbox_dict[b_id][0] = poly.bbox_coord[0]
                if poly.bbox_coord[1] < self.building_bbox_dict[b_id][1]:
                    self.building_bbox_dict[b_id][1] = poly.bbox_coord[1]
                if poly.bbox_coord[2] > self.building_bbox_dict[b_id][2]:
                    self.building_bbox_dict[b_id][2] = poly.bbox_coord[2]
                if poly.bbox_coord[3] > self.building_bbox_dict[b_id][3]:
                    self.building_bbox_dict[b_id][3] = poly.bbox_coord[3]
            id += 1
        self.orig_building_polygon_dict = copy.deepcopy(self.building_polygon_dict)

    def process_gt_rooflines(self, building_id):
        layer = self.file.GetLayerByName("lod22_2d")
        id = 0
        for shape in layer:
            b_id = shape.GetField(building_id)
            poly = WPolygon(id, b_id)
            building_geo = shape.GetGeometryRef()
            if building_geo.GetGeometryType() == ogr.wkbPolygon or \
                    building_geo.GetGeometryName() == "MULTIPOLYGON":
                xmin, xmax, ymin, ymax = building_geo.GetEnvelope()
                poly.bbox_coord = [xmin, ymin, xmax, ymax]  # (lower left, upper right)
                ring = building_geo.GetGeometryRef(0)
                num_points = ring.GetPointCount()
                poly.geo_points = np.empty((num_points, 2), dtype=np.float64)
                for i in range(num_points):
                    x, y, z = ring.GetPoint(i)
                    poly.geo_points[i, 0] = x
                    poly.geo_points[i, 1] = y
            if b_id not in self.building_gt_rooflines_dict.keys():
                self.building_gt_rooflines_dict[b_id] = [poly]
            else:
                self.building_gt_rooflines_dict[b_id].append(poly)
            id += 1

    def preprocess_raw_gpkg(self, args):
        self.srs = self.file.GetLayerByName(args.layer_name).GetSpatialRef()
        self.shapes = self.file.GetLayerByName(args.layer_name)
        if args.merge_connected_building_polygons:
            self.process_to_merge_adjacent_polygons(args.building_id, args.layer_name, args.building_buffer_dis)
        else:
            self.process_raw_shapes(args.building_id, args.layer_name)
        if args.write_gt_building_image_with_rooflines:
            self.process_gt_rooflines(args.building_id)
        del self.shapes
class OrthoPhoto():
    def __init__(self):
        self.image = None
        self.twf_path = None
        self.meta_data = None
        self.pix_size = None
        self.top_left = None
        self.bottom_right = None
        self.srs = None

    def read_image_use_gdal(self, tif_path):
        self.image = gdal.Open(tif_path, gdal.GA_ReadOnly)

    def read_image_meta_from_tif(self):
        extent = self.image.GetGeoTransform()
        self.srs = osr.SpatialReference()
        self.srs.ImportFromWkt(self.image.GetProjection())
        self.pix_size = abs(float(extent[1]))
        self.top_left = (float(extent[0]), float(extent[3]))
        self.bottom_right = (extent[0] + extent[1] * self.image.RasterXSize, extent[3] + extent[5] * self.image.RasterYSize)

    def read_image_meta_from_twf(self, twf_path):
        self.twf_path = twf_path
        tfw_file = open(self.twf_path, "r")
        self.meta_data = tfw_file.read()
        tfw_file.close()
        self.pix_size = abs(float(self.meta_data.split("\n")[3]))
        self.top_left = (float(self.meta_data.split("\n")[4]), float(self.meta_data.split("\n")[5]))
        self.bottom_right = (self.top_left[0] + self.image.RasterXSize * self.pix_size, self.top_left[1] - self.image.RasterYSize * self.pix_size)

class SingleBuildingImage():
    def __init__(self, id, pixel_offset):
        self.id = id
        self.pix_size = None
        self.pixel_offset = pixel_offset
        self.bbox_pix = []
        self.image = None
        self.path = None
        self.xy_start_end = []
        self.srs = None

        self.geo_points = np.zeros((0, 2), dtype=np.float64)
        self.pix_points = np.zeros((0, 2), dtype=np.float64)
        self.edges = np.zeros((0, 2), dtype=np.uint32)
        self.edge_pix_center = []
        self.edge_shape_ind = []

    def sbi_convert_from_geo_to_pix(self, top_left, pix_size, bbox_pix):
        self.pix_points = []
        if len(self.geo_points) != 0:
            for pts in self.geo_points:
                pts_pix = (int((pts[0] - top_left[0]) / pix_size) - bbox_pix[0], int((top_left[1] - pts[1]) / pix_size) - bbox_pix[1])
                self.pix_points.append(pts_pix)
            self.pix_points = np.array(self.pix_points)
        else:
            self.pix_points = np.zeros((0, 2), dtype=np.float64)

    def sbi_convert_from_pix_to_geo(self, top_left, pix_size, bbox_pix):
        self.geo_points = []
        if len(self.pix_points) != 0:
            for pts_pix in self.pix_points:
                pts = (top_left[0] + float((pts_pix[0] + bbox_pix[0]) * pix_size), top_left[1] - float((pts_pix[1] + bbox_pix[1]) * pix_size))
                self.geo_points.append(pts)
            self.geo_points = np.array(self.geo_points)
        else:
            self.geo_points = np.zeros((0, 2), dtype=np.float64)

    def filter_partitions(self, b_poly):
        """
        Keep edges and their corresponding points that are both inside any of the given polygons.
        Return the new edges and points arrays with updated point indices.
        """
        extended_polygons = [ShapeL_polygon(poly.buffered_pix_points) for poly in b_poly]
        kept_points = []  # list to store the points that are kept
        kept_edges = []  # list to store the edges that are kept
        new_point_indices = {}  # dictionary to store new indices of kept points
        next_idx = 0  # counter for new indices of kept points

        # loop through all edges
        for i, (p1_idx, p2_idx) in enumerate(self.edges):
            p1 = self.pix_points[p1_idx]
            p2 = self.pix_points[p2_idx]

            # check if both points are inside any of the polygons
            inside_poly = False
            for poly in extended_polygons:
                if poly.contains(ShapeL_point(p1)) and poly.contains(ShapeL_point(p2)):
                    inside_poly = True
                    break

            # if both points are inside any of the polygons, keep the edge and the points
            if inside_poly:
                if p1_idx not in new_point_indices:
                    new_point_indices[p1_idx] = next_idx
                    kept_points.append(p1)
                    next_idx += 1
                if p2_idx not in new_point_indices:
                    new_point_indices[p2_idx] = next_idx
                    kept_points.append(p2)
                    next_idx += 1
                kept_edges.append([new_point_indices[p1_idx], new_point_indices[p2_idx]])
        # convert kept points and edges to numpy arrays
        if len(kept_points) != 0 and len (kept_edges) != 0:
            self.pix_points = np.array(kept_points, dtype=np.float32)
            self.edges = np.array(kept_edges, dtype=np.int32)
        if self.pix_points.shape[0] == 0 and self.pix_points.shape[0] == 0:
            self.pix_points = np.zeros((0, 2), dtype=np.float64)
            self.edges = np.zeros((0, 2), dtype=np.uint32)

    def collect_centers_and_shapes(self, args, i_building, b_id, building_polys, top_left_used):
        used_shapes = []
        used_polys_pix_center = dict()
        ui = 0
        if args.merge_connected_building_polygons:
            for b_ft_id in building_polys.merged_to_original_map[i_building]:
                used_shapes.extend(building_polys.orig_building_shape_dict[b_ft_id])
                for b_ft_poly in building_polys.orig_building_polygon_dict[b_ft_id]:
                    b_pix_center = b_ft_poly.convert_geo_pt_to_pix_pt(b_ft_poly.geo_center, top_left_used, self.pix_size, self.bbox_pix)
                    if b_ft_id not in used_polys_pix_center.keys():
                        used_polys_pix_center[ui] = [b_pix_center]
                    else:
                        used_polys_pix_center[ui].append(b_pix_center)
                ui += 1
        else:
            used_shapes = building_polys.orig_building_shape_dict[b_id]
            for b_ft_poly in building_polys.orig_building_polygon_dict[b_id]:
                b_pix_center = b_ft_poly.convert_geo_pt_to_pix_pt(b_ft_poly.geo_center, top_left_used, self.pix_size, self.bbox_pix)
                if b_id not in used_polys_pix_center.keys():
                    used_polys_pix_center[ui] = [b_pix_center]
                else:
                    used_polys_pix_center[ui].append(b_pix_center)
                ui += 1

        if len(used_polys_pix_center) > 1:
            for ed in self.edges:
                pt1 = self.pix_points[ed[0]]
                pt2 = self.pix_points[ed[1]]
                midpoint = [(pt1[0] + pt2[0]) / 2, (pt1[1] + pt2[1]) / 2]
                self.edge_pix_center.append(midpoint)

        assert(len(used_shapes) == len(used_polys_pix_center))
        return used_shapes, used_polys_pix_center

    def attach_shape_to_edge(self, used_shapes, used_polys_pix_center):
        if len(used_polys_pix_center) > 1:
            # Build an R-tree index of the points in used_polys_pix_center
            idx = rtree.index.Index()
            for index, points in used_polys_pix_center.items():
                for point in points:
                    idx.insert(index, tuple(point) + tuple(point))

            # Find the nearest point in used_polys_pix_center for each point in edge_pix_center
            self.edge_shape_ind = np.zeros(len(self.edge_pix_center), dtype=int)
            for i, point in enumerate(self.edge_pix_center):
                nearest_index = list(idx.nearest(tuple(point) + tuple(point), 1))[0]
                self.edge_shape_ind[i] = nearest_index
        else:
            self.edge_shape_ind = np.zeros((len(self.edges),), dtype=np.uint32)

    def add_footprint_lines(self, poly):
        # Add new points to the existing array of points
        self.pix_points = np.concatenate((self.pix_points, np.array(poly.pix_points)), axis=0)
        # Get the number of existing points
        num_existing_points = len(self.pix_points) - len(poly.pix_points)
        # Create an array of new edges by adding the index of each new point to the index of the previous point
        new_edges = [[i - 1 + num_existing_points, i + num_existing_points] for i in range(1, len(poly.pix_points))]
        # Add new edges to the existing array of edges
        self.edges = np.concatenate((self.edges, np.array(new_edges, dtype=np.uint32)), axis=0)

    def draw_partitions(self, img):
        if len(self.pix_points) != 0:
            for pts_pix in self.pix_points:
                img = cv2.circle(img, (int(pts_pix[0]), int(pts_pix[1])), 2, (0, 255, 0), 0)
            for ed in self.edges:
                pt1 = self.pix_points[ed[0]]
                pt2 = self.pix_points[ed[1]]
                img = cv2.line(img, (int(pt1[0]), int(pt1[1])), (int(pt2[0]), int(pt2[1])), (0, 0, 255), 1)

    def get_used_images_for_poly(self, img_list, poly_bbox, filename_list):
        corner_img_dict = dict()
        #top_left, top_right, bottom_left, bottom_right
        poly_corners = [(poly_bbox[0], poly_bbox[3]),(poly_bbox[2], poly_bbox[3]), \
                        (poly_bbox[0], poly_bbox[1]),(poly_bbox[2], poly_bbox[1])]
        #get image index and corners index
        img_i = 0
        for img in img_list:
            pi = 0
            for p_corn in poly_corners:
                #check if the corner is in the image
                if p_corn[0] > img.top_left[0] and p_corn[1] < img.top_left[1] \
                        and p_corn[0] < img.bottom_right[0] and p_corn[1] > img.bottom_right[1]:
                    corner_img_dict[pi] = img_i
                pi += 1
            img_i += 1

        #assert len(corner_img_dict.keys()) == 4, f"The input images do not cover whole polygons, missing {4 - len(corner_img_dict.keys())} image(s)"
        if len(corner_img_dict.keys()) != 4:
            return [], [], [], -1
        self.pix_size = img_list[corner_img_dict[0]].pix_size
        self.srs = img_list[corner_img_dict[0]].srs
        # 1
        if corner_img_dict[0] == corner_img_dict[1] \
            and corner_img_dict[0] == corner_img_dict[2] \
            and corner_img_dict[0] == corner_img_dict[3]:
            return img_list[corner_img_dict[0]], img_list[corner_img_dict[0]].top_left, \
                filename_list[corner_img_dict[0]], 0
        # 1
        # 2
        elif corner_img_dict[0] == corner_img_dict[1] \
            and corner_img_dict[2] == corner_img_dict[3] \
            and corner_img_dict[0] != corner_img_dict[2]:
            return np.array([img_list[corner_img_dict[0]], img_list[corner_img_dict[2]]]), \
                img_list[corner_img_dict[0]].top_left, filename_list[corner_img_dict[0]], 1 #vertical
        #1 2
        elif corner_img_dict[0] == corner_img_dict[2] \
            and corner_img_dict[1] == corner_img_dict[3] \
            and corner_img_dict[0] != corner_img_dict[1]:
            return np.array([img_list[corner_img_dict[0]], img_list[corner_img_dict[1]]]),\
                img_list[corner_img_dict[0]].top_left, filename_list[corner_img_dict[0]], 2 #horizontal
        #1 2
        #3 4
        elif corner_img_dict[0] != corner_img_dict[1] \
            and corner_img_dict[0] != corner_img_dict[2] \
            and corner_img_dict[0] != corner_img_dict[3]:
            return np.array([img_list[corner_img_dict[0]], img_list[corner_img_dict[1]], \
                             img_list[corner_img_dict[2]], img_list[corner_img_dict[3]]]), \
                             img_list[corner_img_dict[0]].top_left, filename_list[corner_img_dict[0]], 3

    def set_crop_bbox(self, top_left, poly_bbox):
        self.bbox_pix = []
        self.bbox_pix.append(int((poly_bbox[0] - top_left[0]) / self.pix_size - self.pixel_offset))
        self.bbox_pix.append(int((top_left[1] - poly_bbox[3]) / self.pix_size - self.pixel_offset))
        self.bbox_pix.append(int((poly_bbox[2] - top_left[0]) / self.pix_size + self.pixel_offset))
        self.bbox_pix.append(int((top_left[1] - poly_bbox[1]) / self.pix_size + self.pixel_offset))

    def check_if_offset_exceed_image_boundary(self, top_left, poly_bbox, which_case, img_used):
        offset_limits = 0
        if self.bbox_pix[0] < 0 and abs(self.bbox_pix[0]) > offset_limits:
            offset_limits = abs(self.bbox_pix[0])
        if self.bbox_pix[1] < 0 and abs(self.bbox_pix[1]) > offset_limits:
            offset_limits = abs(self.bbox_pix[1])
        if self.bbox_pix[2] < 0 and abs(self.bbox_pix[2]) > offset_limits:
            offset_limits = abs(self.bbox_pix[2])
        if self.bbox_pix[3] < 0 and abs(self.bbox_pix[3]) > offset_limits:
            offset_limits = abs(self.bbox_pix[3])
        if which_case == 0:
            offset_limits_x = 0
            offset_limits_y = 0
            if self.bbox_pix[2] > img_used.image.RasterXSize:
                offset_limits_x = self.bbox_pix[2] - img_used.image.RasterXSize
            if self.bbox_pix[3] > img_used.image.RasterYSize:
                offset_limits_y = self.bbox_pix[3] - img_used.image.RasterYSize
            if offset_limits_x > offset_limits_y and offset_limits_x > offset_limits:
                offset_limits = offset_limits_x
            elif offset_limits_y > offset_limits_x and offset_limits_y > offset_limits:
                offset_limits = offset_limits_y
        if which_case == 1:
            offset_limits_x = 0
            if self.bbox_pix[2] > img_used[0].image.RasterXSize:
                offset_limits_x = self.bbox_pix[2] - img_used[0].image.RasterXSize
            if offset_limits_x > offset_limits:
                offset_limits = offset_limits_x
        if which_case == 2:
            offset_limits_y = 0
            if self.bbox_pix[3] > img_used[0].image.RasterYSize:
                offset_limits_y = self.bbox_pix[3] - img_used[0].image.RasterYSize
            if offset_limits_y > offset_limits:
                offset_limits = offset_limits_y
        if offset_limits != 0:
            self.pixel_offset = self.pixel_offset - offset_limits
            self.set_crop_bbox(top_left, poly_bbox)

    def crop_image(self, img_used, top_left, poly_bbox, which_case):
        self.set_crop_bbox(top_left, poly_bbox)
        self.check_if_offset_exceed_image_boundary(top_left, poly_bbox, which_case, img_used)
        self.bbox_pix = list(map(int, self.bbox_pix))
        # 1
        if which_case == 0:
            self.image = img_used.image.ReadAsArray(self.bbox_pix[0], self.bbox_pix[1], #x_start, y_start
                                         self.bbox_pix[2] - self.bbox_pix[0], #x size
                                         self.bbox_pix[3] - self.bbox_pix[1]) #y size
            self.image = np.transpose(self.image, (1, 2, 0))
            self.image = cv2.cvtColor(self.image, cv2.COLOR_BGR2RGB)
        # 1
        # 2
        elif which_case == 1:
            img1 = img_used[0].image.ReadAsArray(self.bbox_pix[0], self.bbox_pix[1],
                                      self.bbox_pix[2] - self.bbox_pix[0],
                                      img_used[0].image.RasterYSize - self.bbox_pix[1])
            img2 = img_used[1].image.ReadAsArray(self.bbox_pix[0], 0,
                                      self.bbox_pix[2] - self.bbox_pix[0],
                                      self.bbox_pix[3] - img_used[0].image.RasterYSize)
            img1 = np.transpose(img1, (1, 2, 0))
            img2 = np.transpose(img2, (1, 2, 0))
            img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
            img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
            self.image = np.concatenate((img1, img2), axis=0)
        #1 2
        elif which_case == 2:
            img1 = img_used[0].image.ReadAsArray(self.bbox_pix[0], self.bbox_pix[1],
                                      img_used[0].image.RasterXSize - self.bbox_pix[0],
                                      self.bbox_pix[3] - self.bbox_pix[1])
            img2 = img_used[1].image.ReadAsArray(0, self.bbox_pix[1],
                                      self.bbox_pix[2] - img_used[0].image.RasterXSize,
                                      self.bbox_pix[3] - self.bbox_pix[1])
            img1 = np.transpose(img1, (1, 2, 0))
            img2 = np.transpose(img2, (1, 2, 0))
            img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
            img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
            self.image = np.concatenate((img1, img2), axis=1)
        #1 2
        #3 4
        elif which_case == 3:
            img1 = img_used[0].image.ReadAsArray(self.bbox_pix[0], self.bbox_pix[1],
                                      img_used[0].image.RasterXSize - self.bbox_pix[0],
                                      img_used[0].image.RasterYSize - self.bbox_pix[1])
            img2 = img_used[1].image.ReadAsArray(0, self.bbox_pix[1],
                                      self.bbox_pix[2] - img_used[0].image.RasterXSize,
                                      img_used[1].image.RasterYSize - self.bbox_pix[1])
            img1 = np.transpose(img1, (1, 2, 0))
            img2 = np.transpose(img2, (1, 2, 0))
            img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
            img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
            img1img2 = np.concatenate((img1, img2), axis=1)

            img3 = img_used[0].image.ReadAsArray(self.bbox_pix[0], 0,
                                      img_used[2].image.RasterXSize - self.bbox_pix[0],
                                      self.bbox_pix[3] - img_used[0].image.RasterYSize)
            img4 = img_used[1].image.ReadAsArray(0, 0,
                                      self.bbox_pix[2] - img_used[2].image.RasterXSize,
                                      self.bbox_pix[3] - img_used[1].image.RasterYSize)
            img3 = np.transpose(img3, (1, 2, 0))
            img4 = np.transpose(img4, (1, 2, 0))
            img3 = cv2.cvtColor(img3, cv2.COLOR_BGR2RGB)
            img4 = cv2.cvtColor(img4, cv2.COLOR_BGR2RGB)
            img3img4 = np.concatenate((img3, img4), axis=1)
            self.image = np.concatenate((img1img2, img3img4), axis=0)

    def write_image(self, img_in, path):
        self.path = path
        cv2.imwrite(self.path, img_in)