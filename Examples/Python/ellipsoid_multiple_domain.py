# -*- coding: utf-8 -*-
"""
====================================================================
Full Example Pipeline for Statistical Shape Modeling with ShapeWorks
====================================================================
This example is set to serve as a test case for new ShapeWorks users, and each
step is explained in the shapeworks including the pre-processing, the 
optimization and, the post ShapeWorks visualization.

First import the necessary modules
"""
import os
import glob
import numpy as np
import shapeworks as sw
import OptimizeUtils
import AnalyzeUtils

def Run_Pipeline(args):
    print("\nStep 1. Extract Data\n")
    """
    Step 1: EXTRACT DATA

    We define dataset_name which determines which dataset to download from 
    the portal and the directory to save output from the use case in. 
    """
    print("\nDataset options for running multiple domain use case: \n")
    print("1. ellipsoid_joint_rotation \t 2. ellipsoid_joint_size \t 3. ellipsoid_joint_size_rotation \n")
    print("You can change the dataset name and output directory name to try out this use case with other datasets")


    dataset_name = "ellipsoid_joint_rotation"
    output_directory = "Output/ellipsoid_multiple_domain/"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    # If running a tiny_test, then download subset of the data
    if args.tiny_test:
        sw.data.download_subset(
            args.use_case, dataset_name, output_directory)
        file_list = sorted(glob.glob(output_directory +
                                     dataset_name + "/segmentations/*.nrrd"))[:6]
    # Else download the entire dataset
    else:
        sw.data.download_and_unzip_dataset(dataset_name, output_directory)
        file_list = sorted(glob.glob(output_directory +
                                     dataset_name + "/segmentations/*.nrrd"))

        if args.use_subsample:
            sample_idx = sw.data.sample_images(file_list, int(args.num_subsample),domains_per_shape=2)
            file_list = [file_list[i] for i in sample_idx]

          
    # If skipping grooming, use the pregroomed distance transforms from the portal
    if args.skip_grooming:
        print("Skipping grooming.")
        dt_directory = output_directory + dataset_name + '/groomed/distance_transforms/'
        indices = []
        if args.tiny_test:
            indices = list(range(6))
        elif args.use_subsample:
            indices = sample_idx
        dt_files = sw.data.get_file_list(
            dt_directory, ending=".nrrd", indices=indices)


    # Else groom the segmentations and get distance transforms for optimization
    else:
        print("\nStep 2. Groom - Data Pre-processing\n")
        """
        Step 2: GROOMING 
        The segmentaions are pre-alinged during generation( EllipsoidJointsGenerator) 
        such that they are centered w.r.t to each other. Hence we only perform the 
        following two steps
        The required grooming steps are: 
        1. Isotropic resampling
        2. Select a reference
        3. Rigid alignment
        4. Create smooth signed distance transforms

        For more information on grooming see docs/workflow/groom.md
        http://sciinstitute.github.io/ShapeWorks/workflow/groom.html
        """

        # Create a directory for groomed output
        groom_dir = output_directory + 'groomed/'
        if not os.path.exists(groom_dir):
            os.makedirs(groom_dir)

        """
        First, we need to loop over the shape segmentation files and load the segmentations
        """
        # list of shape segmentations
        shape_seg_list = []
        # list of shape names (shape files prefixes) to be used for saving outputs
        shape_names = []
        domain_ids = []
        for shape_filename in file_list:
            print('Loading: ' + shape_filename)
            # get current shape name
            shape_names.append(shape_filename.split('/')[-1].replace('.nrrd', ''))

            # get domain identifiers
            name = shape_filename.split('/')[-1].replace('.nrrd', '')
            domain_ids.append(name.split(".")[0].split("_")[-1])
            
            # load segmentation
            shape_seg = sw.Image(shape_filename)
            # append to the shape list
            shape_seg_list.append(shape_seg)
        #domain identifiers for all shapes
        domain_ids = np.array(domain_ids)
        #shape index for all shapes in domain 1 
        domain1_indx = list(np.where(domain_ids == 'd1')[0])
        #shape index for all shapes in domain 2
        domain2_indx = list(np.where(domain_ids == 'd2')[0])





        """
        Now we can loop over the segmentations and apply the initial grooming steps to themm
        """
        
        for shape_seg, shape_name in zip(shape_seg_list, shape_names):

            """
            Grooming Step 1: Resample segmentations to have isotropic (uniform) spacing
                - Antialiase the binary segmentation to convert it to a smooth continuous-valued 
                image for interpolation
                - Resample the antialiased image using the same voxel spacing for all dimensions
                - Binarize the resampled images to results in a binary segmentation with an 
                isotropic voxel spacing
            """
            print('Resampling segmentation: ' + shape_name)
            # antialias for 30 iterations
            antialias_iterations = 30
            shape_seg.antialias(antialias_iterations)
            # resample to isotropic spacing using linear interpolation
            iso_spacing = [1, 1, 1]
            shape_seg.resample(iso_spacing, sw.InterpolationType.Linear)
            # make segmetnation binary again
            shape_seg.binarize()

        """
        Grooming Step 2: Select a reference
        This step requires breaking the loop to load all of the segmentations at once so the shape
        closest to the mean can be found and selected as the reference. 
        For multiple domain data, the reference image has to found for the joint as whole. 
        Hence first generate the combined joint 
        """
        resampled_shapes_names = []
        domains_per_shape = 2 
        for shape_seg, shape_name in zip(shape_seg_list,shape_names):
            filename = groom_dir + shape_name.split(".nrrd")[0] + "_resampled.vtk"
            resampled_shapes_names.append(filename)
            shape_seg.toMesh(0.5).write(filename)

        ref_index = sw.find_reference_image_index(resampled_shapes_names,domains_per_shape)
        [os.remove(file) for file in resampled_shapes_names]
        # Make a copy of the reference segmentation
        ref_seg = []
        ref_names = []
        for d in range(domains_per_shape):

            seg = shape_seg_list[domains_per_shape*ref_index + d ].write(groom_dir + 'reference' + str(d).zfill(2)+'.nrrd')
            ref_seg.append(seg.antialias(antialias_iterations))
            ref_names.append(shape_names[domains_per_shape*ref_index + d ])
        print("Reference found: " , ref_names)

        """
        Grooming Step 3: Rigid alignment
        This step rigidly aligns each shape to the selected references. 
        Rigid alignment involves interpolation, hence we need to convert binary segmentations 
        to continuous-valued images again. There are two steps:
            - computing the rigid transformation parameters that would align a segmentation 
            to the reference shape
            - applying the rigid transformation to the segmentation
            - save the aligned images for the next step
        """

        # Set the alignment parameters
        iso_value = 1e-20
        icp_iterations = 200
        for i in range(len(shape_seg_list)):
            ref_domain = ref_seg[i%domains_per_shape]
            ref_name = ref_names[i%domains_per_shape]

            print('Aligning ' + shape_names[i] + ' to ' + ref_name)
            # compute rigid transformation
            shape_seg_list[i].antialias(antialias_iterations)
            rigidTransform = shape_seg_list[i].createTransform(
                ref_domain, sw.TransformType.IterativeClosestPoint, iso_value, icp_iterations)
            # second we apply the computed transformation, note that shape_seg has
            # already been antialiased, so we can directly apply the transformation
            shape_seg.applyTransform(rigidTransform,
                                     ref_domain.origin(),  ref_domain.dims(),
                                     ref_domain.spacing(), ref_domain.coordsys(),
                                     sw.InterpolationType.Linear)
            # then turn antialized-tranformed segmentation to a binary segmentation
            shape_seg_list[i].binarize()

        """
        Grooming Step 2: Converting segmentations to smooth signed distance transforms.
        The computeDT API needs an iso_value that defines the foreground-background interface, to create 
        a smoother interface we first antialiasing the segmentation then compute the distance transform 
        at the zero-level set. We then need to smooth the DT as it will have some remaining aliasing effect 
        of binarization. 
        So the steps are:
            - Antialias 
            - Compute distance transform
            - Apply smoothing
            - Save the distance transform
        """

        # Define distance transform parameters
        iso_value = 0
        sigma = 2
        # Loop over segs and compute smooth DT
        for shape_seg, shape_name in zip(shape_seg_list, shape_names):
            print('Compute DT for segmentation: ' + shape_name)
            shape_seg.antialias(antialias_iterations).computeDT(
                iso_value).gaussianBlur(sigma)
        # Save distance transforms
        dt_files = sw.utils.save_images(groom_dir + 'distance_transforms/', shape_seg_list,
                                        shape_names, extension='nrrd', compressed=False, verbose=True)

    print("\nStep 3. Optimize - Particle Based Optimization\n")
    """
    Step 3: OPTIMIZE - Particle Based Optimization

    Now that we have the distance transform representation of data we create 
    the parameter files for the shapeworks particle optimization routine.
    For more details on the plethora of parameters for shapeworks please refer 
    to docs/workflow/optimze.md
    http://sciinstitute.github.io/ShapeWorks/workflow/optimize.html
    """

    # Make directory to save optimization output
    point_dir = output_directory + 'shape_models/'
    if not os.path.exists(point_dir):
        os.makedirs(point_dir)
    # Create a dictionary for all the parameters required by optimization
    # Create a dictionary for all the parameters required by optimization
    parameter_dictionary = {
        "number_of_particles" : [512,512],
        "use_normals": [0,0],
        "normal_weight": [1.0,1.0],
        "checkpointing_interval" : 200,
        "keep_checkpoints" : 0,
        "iterations_per_split" : 500,
        "optimization_iterations" : 500,
        "starting_regularization" :100,
        "ending_regularization" : 0.5,
        "recompute_regularization_interval" : 2,
        "domains_per_shape" : 2,
        "domain_type" : 'image',
        "relative_weighting" : 1, #10 mesh, # 1 for segmentation images
        "initial_relative_weighting" : 0.1,
        "procrustes_interval" : 0,
        "procrustes_scaling" : 0,
        "save_init_splits" : 0,
        "verbosity" : 3

      }

    if args.tiny_test:
        parameter_dictionary["number_of_particles"] = [32,32]
        parameter_dictionary["optimization_iterations"] = 25
    
    
    [local_point_files, world_point_files] = OptimizeUtils.runShapeWorksOptimize(
        point_dir, dt_files, parameter_dictionary)


    if args.tiny_test:
        print("Done with tiny test")
        exit()

    print("\nStep 4. Analysis - Launch ShapeWorksStudio - sparse correspondence model.\n")
    """
    Step 4: ANALYZE - Shape Analysis and Visualization

    Now we launch studio to analyze the resulting shape model.
    For more information about the analysis step, see docs/workflow/analyze.md
    http://sciinstitute.github.io/ShapeWorks/workflow/analyze.html
    """
    domains_per_shape = 2
    AnalyzeUtils.launchShapeWorksStudio(
        point_dir, dt_files, local_point_files, world_point_files,domains_per_shape)