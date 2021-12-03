import datajoint as dj
import element_data_loader.utils
import pynwb
import warnings

from element_calcium_imaging import scan, imaging

def curated_imaging_to_nwb(curated_imaging_keys,nwbfile):
    """
    Generate one NWBFile object containing all Calcium Imaging Data
        as specified by "curated_imaging_keys"
    Note: specified keys must all be from one session

    :param curated_imaging_keys: entries of imaging.Curation table
    :return: Preexisting NWB file with imaging appended
    """
    ## Validate input
    if isinstance(curated_imaging_keys, dj.expression.QueryExpression):
        curated_imaging_keys = (imaging.Curation & curated_imaging_keys).fetch('KEY')

    assert len(imaging._linking_module.Session & curated_imaging_keys) == 1, \
        f'Multiple sessions error! The specified "curated_imaging_keys"' \
        f' must be from one session only'
    session_key = imaging._linking_module.Session & curated_imaging_keys.fetch('KEY')

    ## Not yet implemented - loop for each provided key
    if len(imaging._linking_module.Curation&curated_imaging_keys) > 1:
        raise NotImplementedError('Multiple imaging keys not yet implemented.\
        Future versions to loop through provided keys.')
    else: curated_imaging_key = curated_imaging_keys

    # ---- Image Specifics ----
    ## Image itself. Chose to load filepath to nwb on first pass
    image_directory = element_data_loader.utils.find_full_path(
        _linking_module.get_imaging_root_data_dir(),
        (session.SessionDirectory&session_key).fetch1('session_dir')
        #OR IS IT: (scan.ScanInfo.ScanFile&session_key).fetch('file_path')
        )
    dimensions=(scan.ScanInfo.Field&curated_imaging_key).fetch(
        'px_width',
        'px_height',
        # 'px_depth' # not yet implemented in element_calcium_imaging
        )
    dimensions=np.stack(dimensions,axis=1) # mult fields, gives [[x,y],[x,y]]
    scale=(scan.ScanInfo.Field&curated_imaging_key).fetch(
        'um_width',
        'um_height',
        # 'um_depth' # not yet implemented in element_calcium_imaging
        )
    if np.isnan(np.sum(scale)): resolution=-1.0 # if any um_ fields blank
    else: resolution=scale/dimensions

    # ---- Equipment Info ----
    ## Not yet implemented in element_lab,
    ##  current implication is fully user specified
    try:
        equipment_info=(lab.Equipement & session_key).fetch1()
        device = nwbfile.create_device(
            name=equipment_info['equipement'],
            description=equipment_info['description'],
            manufacturer=equipment_info['manufacturer']
            )
    except:
        device = None
        warnings.warn('Too many devices or device info not found in lab.Equipement')
    try: # NEEDS WORK - how much info to store in CaImgEquipment, where to put
        CaImg_equip_info = (lab.Equipement.CaImgEquipment & session_key).fetch1()
        optical_channel = OpticalChannel(
            name=CaImg_equip_info['caimg_equip_type'],
            description=CaImg_equip_info['caimg_equip_descrip'],
            emission_lambda=CaImg_equip_info['emission_lambda']
            )
    except: optical_channel=None

    for field in range(int(scan.ScanInfo.fetch('nfields'))):
        imaging_plane = nwbfile.create_imaging_plane(
            name=f'Field {field_info['field_idx']}',
            optical_channel=optical_channel,
            # imaging_rate=30., # usecs_per_line?
            description=field_info['roi'],
            device=device,
            # excitation_lambda=600.,
            # indicator="GFP",
            # location="V1",
            # grid_spacing=[.01, .01],
            # grid_spacing_unit="meters",
            # origin_coords=[1., 2., 3.],
            # origin_coords_unit="meters"
            )

    image_series = pynwb.image.ImageSeries(name='Images',
        # data=None,       # As filepath. Opting not to load in.
        # unit=None,       # Only req if data is specified
        format='external', # As filepath. Opting not to load in
        external_file=image_directory,
        # starting_frame=[0],
        # bits_per_pixel=None,
        dimension=dimensions,
        resolution=resolution,
        # conversion=1.0,
        # timestamps=None,
        # starting_time=None,
        rate=(scan.ScanInfo&session_key).fetch1('fps'),
        comments=(scan.Scan&key).fetch1('scan_notes'),
        description='Structural(?) Images for Scanning',
        # control=None,
        # control_description=None,
        device=device,
        # imaging_plane=imaging_plane
        )
    nwbfile.add_acquisition(image_series)

    # ---- Processed Image ----
    output_directory = element_data_loader.utils.find_full_path(
        _linking_module.get_imaging_root_data_dir(),
        (imaging.Curation&session_key).fetch('file_path')
        )
    ## Motion correction - Rigid
    corrected_rigid = ImageSeries(
        name='corrected',  # this must be named "corrected"
        data=output_directory, # I don't have a different output for rigid and non
        # unit='na',
        format=external,
        # starting_time=0.0,
        # rate=1.0
        )

    xy_translation_rigid = TimeSeries(
        name='xy_translation',
        # data=np.ones((1000, 2)),
        unit='pixels',
        # starting_time=0.0,
        rate=1.0,
    )

    corrected_image_stack_rigid = CorrectedImageStack(
        corrected=corrected_rigid,
        original=image_series,
        xy_translation=xy_translation_rigid,
    )

    ## Motion correction - Nonrigid
    corrected_nonrigid = ImageSeries(
        name='corrected',  # this must be named "corrected"
        # data=output_directory,
        # unit='na',
        # format=external,
        # starting_time=0.0,
        # rate=1.0
        )

    xy_translation_nonrigid = TimeSeries(
        name='xy_translation',
        data=np.ones((1000, 2)),
        unit='pixels',
        # starting_time=0.0,
        rate=1.0,
    )

    corrected_image_stack_nonrigid = CorrectedImageStack(
        corrected=corrected_nonrigid,
        original=image_series,
        xy_translation=xy_translation_nonrigid,
    )
    motion_correction = MotionCorrection(
        corrected_image_stacks=[corrected_image_stack,
                                corrected_image_stack_nonrigid])

    ophys_module = nwbfile.create_processing_module(
        name='ophys',
        description='optical physiology processed data'
    )

    ophys_module.add(motion_correction)

    ## Image segmentation
    img_seg = ImageSegmentation()

    ps = img_seg.create_plane_segmentation(
        name='PlaneSegmentation',
        description='output from segmenting my favorite imaging plane',
        imaging_plane=imaging_plane,
        reference_images=image_series  # optional
    )

    ophys_module.add(img_seg)

    ## ROIs
    for _ in range(30):
        image_mask = np.zeros((100, 100))

        # randomly generate example image masks
        x = np.random.randint(0, 95)
        y = np.random.randint(0, 95)
        image_mask[x:x + 5, y:y + 5] = 1

        # add image mask to plane segmentation
        ps.add_roi(image_mask=image_mask)

    rt_region = ps.create_roi_table_region(
        region=[0, 1],
        description='the first of two ROIs'
    )
    roi_resp_series = RoiResponseSeries(
        name='RoiResponseSeries',
        data=np.ones((50, 2)),  # 50 samples, 2 ROIs
        rois=rt_region,
        unit='lumens',
        rate=30.
    )
    fl = Fluorescence(roi_response_series=roi_resp_series)
    ophys_module.add(fl)