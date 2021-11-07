# -*- coding: utf-8 -*-
"""
====================================================================
Full Example Pipeline for Statistical Shape Modeling with ShapeWorks
====================================================================
The femur data set is comprised of segmented meshes of femurs and corresponding CT
images that are not segmented.
The full images can be carried through every step of grooming.
"""
import os
import glob
import numpy as np
import shapeworks as sw
import OptimizeUtils
import AnalyzeUtils
import json #TODO remove

def Run_Pipeline(args):
    print("\nStep 1. Extract Data\n")
    """
    Step 1: EXTRACT DATA
    We define dataset_name which determines which dataset to download from
    the portal and the directory to save output from the use case in.
    This data is comprised of femur meshes and corresponding unsegmented hip CT scans.
    """
    dataset_name = "femur-v1"
    output_directory = "Output/femur/"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # If running a tiny_test, then download subset of the data
    if args.tiny_test:
        args.use_single_scale = True
        # sw.data.download_subset(args.use_case, dataset_name, output_directory) #TODO uncomment
        mesh_files = sorted(glob.glob(output_directory +
                            dataset_name + "/meshes/*.ply"))[:3]
        image_files = sorted(glob.glob(output_directory +
                            dataset_name + "/images/*.nrrd"))[:3]
        # TODO constraints = 
    # else download the entire dataset
    else:
        sw.data.download_and_unzip_dataset(dataset_name, output_directory)
        mesh_files = sorted(glob.glob(output_directory +
                            dataset_name + "/meshes/*.ply"))
        image_files = sorted(
            glob.glob(output_directory + dataset_name + "/images/*.nrrd"))
        # TODO constraints = 

        # Select data if using subsample
        if args.use_subsample:
            inputMeshes =[sw.Mesh(filename) for filename in mesh_files]
            sample_idx = sw.data.sample_meshes(inputMeshes, int(args.num_subsample))
            mesh_files = [mesh_files[i] for i in sample_idx]
            # TODO constraints =
    
    # TODO remove/fix
    cutting_planes = []
    cutting_plane_counts = []
    constraint_dir = "Output/femur/femur-v1/constraints/"
    for json_file in sorted(os.listdir(constraint_dir)):
        with open(constraint_dir + json_file) as f:
            constraints = json.load(f)
        cutting_planes.append(np.array(constraints["cutting_plane_above_trochanter"]))
        cutting_plane_counts.append(1)
    if args.tiny_test:
        cutting_planes = cutting_planes[:3]
        cutting_plane_counts = cutting_plane_counts[:3]

    # If skipping grooming, use the pregroomed distance transforms from the portal
    if args.skip_grooming:
        print("Skipping grooming.")
        dt_directory = output_directory + dataset_name + '/groomed/distance_transforms/'
        indices = []
        if args.tiny_test:
            indices = [0, 1, 2]
        elif args.use_subsample:
            indices = sample_idx
        dt_files = sw.data.get_file_list(
            dt_directory, ending=".nrrd", indices=indices)
    # Else groom the segmentations and get distance transforms for optimization
    else:
        print("\nStep 2. Groom - Data Pre-processing\n")
        """
        Step 2: GROOMING
        The required grooming steps are:
        1. Apply smoothing and decimation to meshes (Save groomed meshes)
        2. Create a clipped version of meshes
        3. Find reflection transform, apply to clipped (Save reflection transform)
        4. Find COM alignment transform based on clipped (make the center of mass =[0,0,0]), apply to clipped (Save translation)
        5. Select reference mesh from clipped meshes
        6. Find rigid alignment transform using clipped transform w.r.t reference mesh (Save Rigid Transform)
        Option to groom corresponding images (includes applying transforms)
        For more information on grooming see docs/workflow/groom.md
        http://sciinstitute.github.io/ShapeWorks/workflow/groom.html
        """

        # Create a directory for groomed output
        groom_dir = output_directory + 'groomed/'
        if not os.path.exists(groom_dir):
            os.makedirs(groom_dir)

        # Set reference side (arbitrary)
        ref_side = "L" # chosen so reflection happens in tiny test

        """
        To begin grooming, we need to loop over the files and load the meshes and images
        """
        names = []
        mesh_list = []
        clipped_mesh_list = []
        for mesh_filename, cutting_plane in zip(mesh_files, cutting_planes):
            print('\nLoading: ' + mesh_filename)
            # Get shape name
            name = os.path.basename(mesh_filename).replace('.ply', '')
            names.append(name)
            # Get mesh
            mesh = sw.Mesh(mesh_filename)
            mesh_list.append(mesh)

            """
            Grooming Step 1: Smooth and decimate meshes
            """
            print('Smoothing and remeshing: ' + name)
            mesh.smooth(iterations=10).decimate(reduction=0.2).smooth(iterations=10)
            
            """
            Grooming Step 2: Get clipped meshes for alignment
            """
            print('Creating clipped version of: ' + name)
            clipped_mesh = sw.Mesh(mesh_filename).clip(list(cutting_plane)[0], 
                            list(cutting_plane)[1], list(cutting_plane)[2])
            clipped_mesh_list.append(clipped_mesh)

        # Write groomed meshes
        print("\nWriting groomed meshes.")
        mesh_files = sw.utils.save_meshes(groom_dir + 'meshes/', mesh_list,
                            names, extension='vtk', compressed=False, verbose=True)

        # Get alignment transforms
        print("\nFinding initial alignment transforms.")
        reflections = []
        COM_translations = []
        for clipped_mesh, name in zip(clipped_mesh_list, names):
            """
            Grooming Step 2: Get reflection transform - We have left and right femurs, so we reflect the non-reference side
            meshes so that all of the femurs can be aligned.
            """
            print("Finding reflection transform for: " + name)
            if ref_side in name:
                reflection = np.eye(4)
                reflection[0][0] = -1 # Reflect across X
                clipped_mesh.applyTransform(reflection)
            else:
                reflection = np.eye(4) # Identity  
            reflections.append(reflection)

            """
            Grooming Step 3: center of alignment
            Prep for rigid alignment by translating center of mass to [0,0,0]
            """
            print('Finding COM alignment transform ' + name)
            translation = np.eye(4)
            translation[0:3,-1] = -clipped_mesh.centerOfMass()
            clipped_mesh.applyTransform(translation)
            COM_translations.append(translation)

        """
        Grooming Step 3: Select a reference
        This step requires breaking the loop to load all of the shapes at once so the shape
        closest to the mean can be found and selected as the reference.
        """
        print("\nFinding reference.")
        ref_index = sw.find_reference_mesh_index(clipped_mesh_list)
        ref_clipped_mesh = clipped_mesh_list[ref_index].write(groom_dir + 'reference.vtk')
        ref_name = names[ref_index]
        print("Reference found: " + ref_name + "\n")

        """
        Grooming Step 4: Rigid alignment
        This step rigidly aligns each shape to the selected references.
        """
        rigid_transforms = []
        for clipped_mesh, name in zip(clipped_mesh_list, names):
            print('Finding alignment transform from ' + name + ' to ' + ref_name)
            rigid_transform = clipped_mesh.createTransform(ref_clipped_mesh, sw.Mesh.AlignmentType.Rigid)
            clipped_mesh.applyTransform(rigid_transform)
            rigid_transforms.append(rigid_transform)
    
        """
        Groom images
        """
        if args.groom_images:

            # Load corresponding images
            print("\nLoading images.")
            image_list = []
            for name in names:
                # Get corresponding image path
                prefix = name.split("_")[0]
                for index in range(len(image_files)):
                    if prefix in image_files[index]:
                        corresponding_image_file = image_files[index]
                        break
                print('Loading image: ' + name)
                image = sw.Image(corresponding_image_file)
                image_list.append(image)

            # Apply transforms to images
            bounding_box = sw.MeshUtils.boundingBox(clipped_mesh_list)
            print(bounding_box)
            for image, name, reflection, COM_translation, rigid_transform in zip(image_list, names, reflections, COM_translations, rigid_transforms):
                print("\nGrooming image: " + name)
                print("Reflecting image: " + name)
                image.applyTransform(reflection)
                print("Translating image: " + name)
                image.setOrigin(image.origin() + COM_translation[0:3,-1])
                print("Aligning image: " + name)
                image.applyTransform(rigid_transform)
                print('Cropping image: ' + name)
                image.crop(bounding_box)

            # Write images
            print("\nWriting groomed images.")
            image_files = sw.utils.save_images(groom_dir + 'images/', image_list,
                            names, extension='nrrd', compressed=True, verbose=True)

        # TODO remove
        print("\nTemporarily writing clipped meshes.")
        mesh_files = sw.utils.save_meshes(groom_dir + 'clipped_meshes/', clipped_mesh_list,
                            names, extension='vtk', compressed=False, verbose=True)

    print("\nStep 3. Optimize - Particle Based Optimization\n")
    # TODO pass transforms
    """
    Step 3: OPTIMIZE - Particle Based Optimization
    Now that we have the distance transform representation of data we create 
    the parameter files for the shapeworks particle optimization routine.
    For more details on the plethora of parameters for shapeworks please refer 
    to docs/workflow/optimze.md
    http://sciinstitute.github.io/ShapeWorks/workflow/optimize.html
    """

    # Make directory to save optimization output
    point_dir = output_directory + 'shape_models/' + args.option_set
    if not os.path.exists(point_dir):
        os.makedirs(point_dir)
    # Create a dictionary for all the parameters required by optimization
    parameter_dictionary = {
        "number_of_particles": 512,
        "use_normals": 0,
        "normal_weight": 10.0,
        "checkpointing_interval": 200,
        "keep_checkpoints": 0,
        "iterations_per_split": 1000,
        "optimization_iterations": 500,
        "starting_regularization": 100,
        "ending_regularization": 0.1,
        "recompute_regularization_interval": 2,
        "domains_per_shape" : 1,
        "domain_type" : 'mesh',
        "relative_weighting" : 10,
        "initial_relative_weighting" : 0.01,
        "procrustes_interval" : 1,
        "procrustes_scaling" : 1,
        "save_init_splits" : 1,
        "verbosity" : 0,
        "use_statistics_in_init" : 0
        # "cutting_plane_counts": cutting_plane_counts, #TODO uncomment when grooming transforms are passed instead of applied
        # "cutting_planes": cutting_planes
    }
     # If running a tiny test, reduce some parameters
    if args.tiny_test:
        parameter_dictionary["number_of_particles"] = 32
        parameter_dictionary["optimization_iterations"] = 25
        parameter_dictionary["iterations_per_split"] = 25
    # Run multiscale optimization unless single scale is specified
    if not args.use_single_scale:
        parameter_dictionary["use_shape_statistics_after"] = 64

    # Execute the optimization function on distance transforms
    [local_point_files, world_point_files] = OptimizeUtils.runShapeWorksOptimize(
        point_dir, mesh_files, parameter_dictionary)

    # If tiny test or verify, check results and exit
    AnalyzeUtils.check_results(args, world_point_files) #TODO uncomment

    print("\nStep 4. Analysis - Launch ShapeWorksStudio - sparse correspondence model.\n")
    """
    Step 4: ANALYZE - Shape Analysis and Visualization
    Now we launch studio to analyze the resulting shape model.
    For more information about the analysis step, see docs/workflow/analyze.md
    http://sciinstitute.github.io/ShapeWorks/workflow/analyze.html
    """
    # Prepare analysis XML
    analyze_xml = point_dir + "/femur_analyze.xml"
    AnalyzeUtils.create_analyze_xml(analyze_xml, mesh_files, local_point_files, world_point_files)
    AnalyzeUtils.launch_shapeworks_studio(analyze_xml)