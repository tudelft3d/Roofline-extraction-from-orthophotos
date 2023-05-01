import copy
import sys
import os.path
import glob
import argparse
import class_definition as cd
import random
import copy
import time
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from joblib import parallel_backend
#from memory_profiler import profile

#modify it to your build path
sys.path.append("./kinetic_partition/build")
import libkinetic_partition
# Define the function to run in parallel

#parse the arguments #
parser = argparse.ArgumentParser(description='Generate rooflines from orthophoto for 3D BAG building model reconstruction.')

# IO configuration
parser.add_argument('--ROOT_PATH', default='/media/geo3d/data/rooflines_data/data')
parser.add_argument('--tif_with_tfw', default=False, type=bool,
                    help='True for with tfw; False without tfw.')
parser.add_argument('--layer_name', default='pand', help='Used layer name in the gis data.') #lod22_2d
parser.add_argument('--building_id', default='fid', help='Used field id for representing individual buildings in the gis data.')
parser.add_argument('--write_building_image', default=True, type=bool,
                    help='Write cropped building image.')
parser.add_argument('--write_building_image_with_rooflines', default=True, type=bool,
                    help='Write image with projected rooflines.')
parser.add_argument('--write_gt_building_image_with_footprints', default=2, type=int,
                    help='Write image with projected ground truth (3D BAG pand) footprints:'
                         '0: not write'
                         '1: write original footprint'
                         '2: write buffered footprint or merged footprint or both.')
parser.add_argument('--write_gt_building_image_with_rooflines', default=True, type=bool,
                    help='Write image with projected ground truth (3D BAG lod22_2d) rooflines.')
parser.add_argument('--multi_processing_cores', default=True, type=int,
                    help='Write image with projected ground truth (3D BAG lod22_2d) rooflines.')

# Processing parameters for building image and polygons
parser.add_argument('--apply_polygon_buffer_filter', default=True, type=bool,
                    help='Use the extended footprint to crop the partition results.')
parser.add_argument('--merge_connected_building_polygons', default=True, type=bool,
                    help='Merge connected building polygons to increase the building completeness.')
parser.add_argument('--add_footprint_to_rooflines', default=True, type=bool,
                    help='Add footprint lines to rooflines for building completeness.')
parser.add_argument('--polygon_buffer', default=60, type=float,
                    help='Extend the polygon with buffer (pixels) to crop the partition results.')
parser.add_argument('--pixel_offset', default=100, type=int,
                    help='Pixel offset for crop the building image')
parser.add_argument('--building_buffer_dis', default=1.0, type=float,
                    help='The minimum distance (> 0) between adjacent buildings for merging.')
parser.add_argument('--show_progress', default=True, type=bool, help='Show progress bar.')
parser.add_argument('--parallel_processing', default=True, type=bool,
                    help='Allows for parallel processing.')
parser.add_argument('--num_cpu_cores_to_use', default=10, type=int, help='Number of CPU cores to use.')
parser.add_argument('--parallel_verbose', default=1, type=int, help='Controls the level of progress reporting.')

# kinetic partition parameters
parser.add_argument('--lsd_scale', default=0.8, type=float,
                    help='Scale the image by Gaussian filter to scale for line detection.')
parser.add_argument('--num_intersection', default=1, type=int,
                    help='Intersecting times of kinetic partition.')
parser.add_argument('--enable_regularise', default=True, type=bool,
                    help='Enable regularization for the detected lines.')
parser.add_argument('--verbose', default=False, type=bool,
                    help='Print each steps for debugging.')
args = parser.parse_args()
#@profile
def read_ortho_photos(ortho_images_in, input_img_dir, args):
    file_name_list = []
    ortho_image_list = []
    n_ortho_images = len(ortho_images_in)
    i_file = 0
    for image_in in ortho_images_in:
        file_name = os.path.splitext(os.path.basename(image_in))[0]
        file_name_suffix = file_name.split("_")
        file_name_used = file_name_suffix[len(file_name_suffix) - 2] + "_" + file_name_suffix[len(file_name_suffix) - 1]
        file_name_list.append(file_name_used)
        meta_in = input_img_dir + file_name + ".tfw"

        i_file = i_file + 1
        #print("Read " + str(i_file) + " / " + str(n_ortho_images) + " ---> " + file_name)

        ortho_image = cd.OrthoPhoto()
        ortho_image.read_image_use_gdal(image_in)
        if args.tif_with_tfw:
            ortho_image.read_image_meta_from_twf(meta_in)
        else:
            ortho_image.read_image_meta_from_tif()
        ortho_image_list.append(ortho_image)
    return ortho_image_list, file_name_list

#@profile
def building_rooflines_generation_pipeline(args, polygon_files_in, i_file, output_rooflines_dir, ortho_image_list,
                                           file_name_list, output_img_dir, output_img_gt_footprint_dir, output_img_rooflines_dir,
                                           output_img_gt_rooflines_dir, n_poly_files):
    poly_in = polygon_files_in[i_file]
    poly_name = os.path.splitext(os.path.basename(poly_in))[0]
    if not args.parallel_processing and not args.show_progress:
        print("Read " + str(i_file) + " / " + str(n_poly_files) + " ---> " + poly_name)
    building_polys = cd.WPolygons()
    building_polys.read_poly(poly_in, args)

    # initialize the polygons to write if there is predict data
    rooflines_out = cd.WPolygons()
    rooflines_out.initialize_gpkg_for_write(output_rooflines_dir + poly_name + "_rooflines", building_polys)

    i_building = 0
    n_building_files = len(building_polys.building_polygon_dict.items())
    # Loop all building polygons to crop image
    for b_id, b_poly in building_polys.building_polygon_dict.items():
        if not args.parallel_processing and not args.show_progress:
            print("Generating building " + str(i_building) + " / " + str(n_building_files) + " ---> " + str(b_id))
        # Extract single builing image
        sing_b = cd.SingleBuildingImage(b_id, args.pixel_offset)
        img_used, top_left_used, image_name_used, which_case = sing_b.get_used_images_for_poly(ortho_image_list,
                                                               building_polys.building_bbox_dict[b_id], file_name_list)
        if isinstance(img_used, list) and img_used == []:
            continue
        sing_b.crop_image(img_used, top_left_used, building_polys.building_bbox_dict[b_id], which_case)
        sing_b.write_image(sing_b.image, output_img_dir + image_name_used + "_" + str(b_id) + ".jpg")

        # write 2d footprints
        sing_b_img_gt_footprint = copy.deepcopy(sing_b.image)
        for poly in b_poly:
            poly.convert_from_geo_to_pix(top_left_used, sing_b.pix_size, sing_b.bbox_pix)
            if args.apply_polygon_buffer_filter:
                poly.buffering_polygon(args)
            if args.write_gt_building_image_with_footprints > 1 or \
                    (args.write_gt_building_image_with_footprints == 1 and \
                     not args.apply_polygon_buffer_filter and not args.merge_connected_building_polygons):
                poly.draw_polygon(sing_b_img_gt_footprint, args.write_gt_building_image_with_footprints)

        # perform partition to extract rooflines
        img_path_read = output_img_dir + image_name_used + "_" + str(b_id) + ".jpg"
        sing_b.pix_points, sing_b.edges = libkinetic_partition.partition_image(img_path_read, args.lsd_scale,
                                                                               args.num_intersection,
                                                                               args.enable_regularise, args.verbose)

        # Fitler out some points and edges based on the input 2d footprints
        if args.apply_polygon_buffer_filter:
            sing_b.filter_partitions(b_poly)

        # Add points and edges from input 2d footprints
        if args.add_footprint_to_rooflines:
            if args.merge_connected_building_polygons:
                for b_ft_id in building_polys.merged_to_original_map[i_building]:
                    for b_ft_poly in building_polys.orig_building_polygon_dict[b_ft_id]:
                        b_ft_poly.convert_from_geo_to_pix(top_left_used, sing_b.pix_size, sing_b.bbox_pix)
                        sing_b.add_footprint_lines(b_ft_poly)
                        if args.write_gt_building_image_with_footprints == 1:
                            b_ft_poly.draw_polygon(sing_b_img_gt_footprint,
                                                   args.write_gt_building_image_with_footprints)
            else:
                for poly in b_poly:
                    sing_b.add_footprint_lines(poly)

        if not args.write_building_image:
            os.remove(img_path_read)

        # Assign each edge a shape for attribute parsing
        used_shapes, used_polys_pix_center = sing_b.collect_centers_and_shapes(args, i_building, b_id, building_polys, top_left_used)
        sing_b.attach_shape_to_edge(used_shapes, used_polys_pix_center)

        # write footprint
        if args.write_gt_building_image_with_footprints > 0:
            sing_b.write_image(sing_b_img_gt_footprint,
                               output_img_gt_footprint_dir + image_name_used + "_" + str(b_id) + "_footprint.jpg")
        del sing_b_img_gt_footprint

        # write results to images
        sing_b_img_poly = copy.deepcopy(sing_b.image)
        sing_b.draw_partitions(sing_b_img_poly)
        sing_b.sbi_convert_from_pix_to_geo(top_left_used, sing_b.pix_size, sing_b.bbox_pix)
        if args.write_building_image_with_rooflines:
            sing_b.write_image(sing_b_img_poly,
                               output_img_rooflines_dir + image_name_used + "_" + str(b_id) + "_rooflines.jpg")
        del sing_b_img_poly

        # write ground truth rooflines if set to true
        sing_b_img_poly_gt = copy.deepcopy(sing_b.image)
        if args.write_gt_building_image_with_rooflines:
            if args.merge_connected_building_polygons:
                for b_gt_id in building_polys.merged_to_original_map[i_building]:
                    if b_gt_id in building_polys.building_gt_rooflines_dict:
                        for b_gt_poly in building_polys.building_gt_rooflines_dict[b_gt_id]:
                            b_gt_poly.convert_from_geo_to_pix(top_left_used, sing_b.pix_size, sing_b.bbox_pix)
                            b_gt_poly.draw_polygon(sing_b_img_poly_gt, 1)
            else:
                if b_id in building_polys.building_gt_rooflines_dict:
                    for b_gt_poly in building_polys.building_gt_rooflines_dict[b_id]:
                        b_gt_poly.convert_from_geo_to_pix(top_left_used, sing_b.pix_size, sing_b.bbox_pix)
                        b_gt_poly.draw_polygon(sing_b_img_poly_gt, 1)
            sing_b.write_image(sing_b_img_poly_gt,
                               output_img_gt_rooflines_dir + image_name_used + "_" + str(b_id) + "_rooflines_gt.jpg")
        del sing_b_img_poly_gt

        # write lines to gpkg
        rooflines_out.add_partition_lines(b_id, sing_b, used_shapes)

        i_building = i_building + 1
        del sing_b.image
        del sing_b.pix_points
        del sing_b.edges
        del sing_b
        del img_used
        del top_left_used
        del image_name_used
        del which_case
        del img_path_read
    rooflines_out.close_file(args)
    building_polys.close_file(args)
    del building_polys
    if args.parallel_processing:
        del ortho_image_list
        del file_name_list

if __name__ == '__main__':
    print(args.ROOT_PATH)
    if args.parallel_processing:
        print("Using " + str(args.num_cpu_cores_to_use) + " cpu cores for parallel processing.")
    root = args.ROOT_PATH + '/'
    input_polygons_dir = root + "poly/" #2D BAG building footprints
    input_img_dir = root + "img/"

    if not os.path.isdir(input_polygons_dir):
        raise ValueError("%s does not exist" % input_polygons_dir)

    if not os.path.isdir(input_img_dir):
        raise ValueError("%s does not exist" % input_img_dir)

    if not os.path.isdir(root + "rgb"):
        os.mkdir(root + "rgb")
    if not os.path.isdir(root + "rooflines"):
        os.mkdir(root + "rooflines")
    if not os.path.isdir(root + "rgb_rooflines"):
        os.mkdir(root + "rgb_rooflines")

    if not os.path.isdir(root + "rgb_gt_footprint") and args.write_gt_building_image_with_footprints > 0:
        os.mkdir(root + "rgb_gt_footprint")
    if not os.path.isdir(root + "rgb_gt_rooflines") and args.write_gt_building_image_with_rooflines > 0:
        os.mkdir(root + "rgb_gt_rooflines")

    output_img_dir = root + "rgb/"
    output_rooflines_dir = root + "rooflines/"
    output_img_rooflines_dir = root + "rgb_rooflines/"

    output_img_gt_footprint_dir = root + "rgb_gt_footprint/"
    output_img_gt_rooflines_dir = root + "rgb_gt_rooflines/"

    #Generate ground-truth from 3D BAG polyfile and orthophotos
    polygon_files_in = glob.glob(input_polygons_dir + "*.gpkg")
    ortho_images_in = glob.glob(input_img_dir + "*.tif")
    ortho_metas_in = glob.glob(input_img_dir + "*.tfw")

    if (len(polygon_files_in) == 0):
        raise ValueError('%s is empty' % input_polygons_dir)
    if (len(ortho_images_in) == 0):
        raise ValueError('%s is empty' % input_img_dir)
    if args.tif_with_tfw:
        if (len(ortho_images_in) != len(ortho_metas_in)):
            raise ValueError('The input images *.tif and meta data *.tfw in %s are not all match' %(input_img_dir))

    #read polyfile
    start_time = time.time()
    n_poly_files = len(polygon_files_in)
    if args.parallel_processing:
        # Run the function in parallel for each file
        # Define a function that calls building_rooflines_generation_pipeline
        def process_file(i_file):
            ortho_image_list, file_name_list = read_ortho_photos(ortho_images_in, input_img_dir, args)
            return building_rooflines_generation_pipeline(args, polygon_files_in, i_file, output_rooflines_dir,
                                                             ortho_image_list, file_name_list, output_img_dir,
                                                             output_img_gt_footprint_dir, output_img_rooflines_dir,
                                                             output_img_gt_rooflines_dir, n_poly_files)
        if args.show_progress:
            with parallel_backend('multiprocessing'):
                results = Parallel(n_jobs=args.num_cpu_cores_to_use, verbose=0)(
                    delayed(process_file)(i_file) for i_file in tqdm(range(n_poly_files), desc="Processing Files", miniters=1))
        else:
            with parallel_backend('multiprocessing'):
                results = Parallel(n_jobs=args.num_cpu_cores_to_use, verbose=args.parallel_verbose)(
                    delayed(process_file)(i_file) for i_file in range(n_poly_files))
    else:
        # read ortho photos
        print("Read all ortho photos.")
        ortho_image_list, file_name_list = read_ortho_photos(ortho_images_in, input_img_dir, args)

        if args.show_progress:
            for i_file in tqdm(range(n_poly_files)):
                building_rooflines_generation_pipeline(args, polygon_files_in, i_file, output_rooflines_dir, ortho_image_list,
                                                          file_name_list, output_img_dir, output_img_gt_footprint_dir,
                                                          output_img_rooflines_dir, output_img_gt_rooflines_dir, n_poly_files)
        else:
            for i_file in range(n_poly_files):
                building_rooflines_generation_pipeline(args, polygon_files_in, i_file, output_rooflines_dir, ortho_image_list,
                                                          file_name_list, output_img_dir, output_img_gt_footprint_dir,
                                                          output_img_rooflines_dir, output_img_gt_rooflines_dir, n_poly_files)

    end_time = time.time()
    processing_time = end_time - start_time
    print(f"Processing time: {processing_time} seconds")