# Laser Information
# Velodyne HDL-64E S2
# Wavelength: 905 nm
# Pulselength: 5 nanoseconds
# Angular resolution: 0.09 degree
# Distance accuracy (sigma): 2cm
# Range: 50m (~0.1 reflectivity), 120m (~0.8 reflectivity)

# Laser Information for Velodyne VLP-16:
# Wavelength: 903 nm
# Pulse length: 6 nanoseconds
# Angular resolution (horizontal): 0.1 - 0.4 deg
# Angular resolution (vertical): 2 deg
# Distance accuracy (typical, sigma): 3cm
# Range: up to 100m
# Rotation: between 300 and 1200 RPM with 60 RPM increments
# the angular resolution in the azimuth depends on the RPM (see manual)

import math
import sys
import os
import struct 
import ctypes
import time 
import random 
import bpy
from mathutils import Vector, Euler, Matrix

from blensor import evd
from blensor import mesh_utils

import blensor
import numpy

BLENSOR_VELODYNE_HDL64E2 = "hdl64e2"
BLENSOR_VELODYNE_HDL32E = "hdl32e"
BLENSOR_VELODYNE_VLP16 = "vlp16"



parameters = {"angle_resolution":0.1728, "rotation_speed":10,"max_dist":120,"noise_mu":0.0,"noise_sigma":0.01,
              "start_angle":0,"end_angle":360, "distance_bias_noise_mu": 0, "distance_bias_noise_sigma": 0.078,
              "reflectivity_distance":50,"reflectivity_limit":0.1,"reflectivity_slope":0.01,
              "noise_types": [("gaussian", "Gaussian", "Gaussian distribution (mu/simga)"),("laplace","Laplace","Laplace distribution (sigma=b)")],
              "models": [(BLENSOR_VELODYNE_HDL64E2, "HDL-64E2", "HDL-64E2"), (BLENSOR_VELODYNE_HDL32E, "HDL-32E", "HDL-32E"), (BLENSOR_VELODYNE_VLP16, "VLP-16", "VLP-16")],
              "output_laser_id_as_color": False}


vlp16_parameters = {
    "angle_resolution": 0.2,
    "rotation_speed": 10, # in Hz, equivalent to 600 RPM
    "max_dist": 100,
    "noise_mu": 0.0,
    "noise_sigma": 0.01,
    "start_angle": 0,
    "end_angle": 360,
    "distance_bias_noise_mu": 0,
    "distance_bias_noise_sigma": 0.078,
    "reflectivity_distance": 50,
    "reflectivity_limit": 0.1,
    "reflectivity_slope": 0.01,
    "noise_types": [("gaussian", "Gaussian", "Gaussian distribution (mu/simga)"),
                    ("laplace", "Laplace", "Laplace distribution (sigma=b)")],
    "models": [(BLENSOR_VELODYNE_HDL64E2, "HDL-64E2", "HDL-64E2"), (BLENSOR_VELODYNE_HDL32E, "HDL-32E", "HDL-32E"),
               (BLENSOR_VELODYNE_VLP16, "VLP-16", "VLP-16")]
}

def addProperties(cType):
    global parameters
    cType.velodyne_angle_resolution = bpy.props.FloatProperty( name = "Scan resolution", default = parameters["angle_resolution"], min = 0.01, max = 3.141, description = "How far(angle) two scan lines are apart" )
    cType.velodyne_rotation_speed = bpy.props.FloatProperty( name = "Rotation speed", default = parameters["rotation_speed"], min = 1, max = 100, description = "Rotation speed in Hertz" )
    cType.velodyne_max_dist = bpy.props.FloatProperty( name = "Scan distance", default = parameters["max_dist"], min = 0, max = 1000, description = "How far the laser can see" )
    cType.velodyne_noise_mu = bpy.props.FloatProperty( name = "Noise mu", default = parameters["noise_mu"], description = "The center of the gaussian noise" )
    cType.velodyne_noise_sigma = bpy.props.FloatProperty( name = "Noise sigma", default = parameters["noise_sigma"], description = "The sigma of the gaussian noise" )
    cType.velodyne_db_noise_mu = bpy.props.FloatProperty( name = "DB Noise mu", default = parameters["distance_bias_noise_mu"], description = "The center of the gaussian noise" )
    cType.velodyne_db_noise_sigma = bpy.props.FloatProperty( name = "DB Noise sigma", default = parameters["distance_bias_noise_sigma"], description = "The sigma of the gaussian noise" )
    cType.velodyne_start_angle = bpy.props.FloatProperty( name = "Start angle", default = parameters["start_angle"], description = "The angle at which the scan is started" )
    cType.velodyne_end_angle = bpy.props.FloatProperty( name = "End angle", default = parameters["end_angle"], description = "The angle at which the scan is stopped" )
 
    cType.velodyne_ref_dist = bpy.props.FloatProperty( name = "Reflectivity Distance", default = parameters["reflectivity_distance"], description = "Objects closer than reflectivity distance are independent of their reflectivity" )
    cType.velodyne_ref_limit = bpy.props.FloatProperty( name = "Reflectivity Limit", default = parameters["reflectivity_limit"], description = "Minimum reflectivity for objects at the reflectivity distance" )
    cType.velodyne_ref_slope = bpy.props.FloatProperty( name = "Reflectivity Slope", default = parameters["reflectivity_slope"], description = "Slope of the reflectivity limit curve" )
 
    cType.velodyne_noise_type = bpy.props.EnumProperty( items= parameters["noise_types"], name = "Noise distribution", description = "Which noise model to use for the distance bias" )
    cType.velodyne_model = bpy.props.EnumProperty( items= parameters["models"], name = "Model", description = "Velodyne Model" )
    cType.velodyne_output_laser_id_as_color = bpy.props.BoolProperty(default=parameters["output_laser_id_as_color"],
                                                                     name="Output laser id as color",
                                                                     description="If enabled, the laser ids will be returned as the color of a sample")
    # TODO add the new parameters here as well

 


def deg2rad(deg):
    return deg*math.pi/180.0

def rad2deg(rad):
    return rad*180.0/math.pi

def tuples_to_list(tuples):
    l = []
    for t in tuples:
        l.extend(t)
    return l


laser_angles =[-7.1143909000 ,-6.8259001000 ,0.3328709900 ,0.6607859700 ,
-6.4908152000 ,-6.0973902000 ,-8.5282297000 ,-8.1613369000 ,-5.8425112000 ,
-5.4713659000 ,-7.8512449000 ,-7.5291781000 ,-3.0490510000 ,-2.7686839000 ,
-5.0532770000 ,-4.7975101000 ,-2.3829660000 ,-2.1140020000 ,-4.4367881000 ,
-4.0757122000 ,-1.7279370000 ,-1.4470620000 ,-3.7842851000 ,-3.4226439000 ,
1.0237820000 ,1.3398750000 ,-0.9904980100 ,-0.6977589700 ,1.6909920000 ,
1.9717960000 ,-0.3464250000 ,-0.0419129990 ,-22.6174770000 ,-22.2171250000 ,
-11.3515270000 ,-10.7762110000 ,-21.7437740000 ,-21.1960530000 ,-24.8932860000 ,
-24.2637140000 ,-20.6647490000 ,-20.0160450000 ,-23.7457070000 ,-23.0651020000 ,
-16.4140130000 ,-16.0144540000 ,-19.4244710000 ,-19.1009100000 ,-15.4937170000 ,
-15.0140580000 ,-18.6500020000 ,-18.1543940000 ,-14.4663670000 ,-13.8276510000 ,
-17.5921270000 ,-16.9942110000 ,-10.3347680000 ,-9.8352394000 ,-13.2298120000 ,
-12.8963990000 ,-9.3798056000 ,-8.8888798000 ,-12.3722690000 ,-11.9693750000]

laser_angles_32 = [-30.67, -9.33, -29.33, -8.00, -28.00, -6.66, -26.66,
                   -5.33, -25.33, -4.00, -24.00, -2.67, -22.67, -1.33,
                   -21.33, 0.00, -20.00, 1.33, -18.67, 2.67, -17.33,
                   4.00, -16.00, 5.33, -14.67, 6.67, -13.33, 8.00,
                   -12.00, 9.33, -10.67, 10.67]

# laser ID is index in list
laser_angles_vlp16 = [-15, 1, -13, 3, -11, 5, -9, 7, -7, 9, -5, 11, -3, 13, -1, 15]


# in mm, see VLP-16 manual, these should be added to the z-axis
vertical_corrections_vlp16 = [11.2, -0.7, 9.7, -2.2, 8.1, -3.7, 6.6, -5.1, 5.1, -6.6, 3.7, -8.1, 2.2, -9.7, 0.7, -11.2]


# The laser noise is initialized with a fixed randomized array to increase
# reproducibility. If the noise should be randomize, call 
# randomize_distance_bias
laser_noise =  [0.023188431056485468, 0.018160539830319688, 
               -0.082857233607375583, -0.064524698320918547, 
               -0.093271246114693618, 0.043527643100361682,
                0.02964170651252759, 0.022130151884228375, 
                0.057142456134798118, -0.0022902380317122851, 
                0.0067235936679492184, -0.036271973173879445, 
                0.017223035592365758, -0.077567037832670549, 
               -0.045212530018321921, -0.071407355488847954, 
               -0.044010545427414234, -0.13100196824243787, 
                0.06562213285464942, 0.072876545417746588, 
                0.02948486101287804, -0.0549980634164711, 
                0.0017215074797752329, 0.011024793340907989, 
               -0.031627028532767984, -0.0015734962901194415, 
               -0.036013321659115423, 0.10758083132777119, 
               -0.08155406184051478, 0.042643613869132221, 
               -0.0084507159317202072, -0.13509680354593059, 
               -0.011626407663142539, -0.016650248443238872, 
               -0.089390919093335297, 0.058944202109090765, 
               -0.0017922326055001933, 0.17981608947109032, 
               -0.08508532554265108, 0.073143736372048768, 
               -0.048115410987441362, -0.042704861239692658, 
                0.035086496817459518, 0.064943607611183743, 
                0.010321640735410247, -0.088027285788537329, 
                0.064927912617182712, -0.0063073544055394313, 
                0.0094236036091259051, -0.012188267349634907, 
               -0.030570059320519608, -0.022883795755219431, 
               -0.01536736072143652, 0.091570972465821771, 
               -0.09083766812058082, 0.14380982060713934, 
                0.01373428104250375, -0.014464880110123654, 
               -0.031761209027541266, -0.01571113598069827, 
               -0.14107381715154735, -0.064936750764770235, 
                0.034911820770082327, 0.065682492298063416]


# polynomial coefficients for noise stddevs from 1m to 5m distance
# rows are lasers, index being laser id
wall_noise_polynomials = list(map(numpy.poly1d, [[1.149959777101233717e-04, -1.447756815982880786e-03, 6.510018018593014509e-03,
                           -1.223405369518627307e-02, 1.650758676374098829e-02],
                          [1.020163480950387355e-04, -1.347023637972283002e-03, 6.462602470486273655e-03,
                           -1.305210838199708477e-02, 1.784158289945754483e-02],
                          [1.146020230467176357e-04, -1.393072554819377697e-03, 5.792563689378243431e-03,
                           -9.648997906561112073e-03, 1.417917061964388097e-02],
                          [8.650362566993354159e-05, -1.135689903779567951e-03, 5.464422232895918725e-03,
                           -1.136699839596469301e-02, 1.736016550102446182e-02],
                          [6.389678743342807265e-05, -6.941993600298232816e-04, 2.840338958273228107e-03,
                           -5.360543336878565573e-03, 1.221417697163846086e-02],
                          [-2.064130027892076505e-05, 1.114051046576767660e-04, 3.967506880809539439e-04,
                           -2.847855223363543638e-03, 1.255669917294362052e-02],
                          [1.307615087548785748e-04, -1.384753774386732975e-03, 5.311792326053499380e-03,
                           -8.893292093065656434e-03, 1.425053945508668145e-02],
                          [-2.928430250961752750e-05, 2.854128538385984279e-04, -4.614062349368840002e-04,
                           -1.716008808146514972e-03, 1.338855253777419992e-02],
                          [-1.706775798030039947e-04, 1.864701880524312827e-03, -7.000606371725321894e-03,
                           1.059679020136581468e-02, 4.753019794913990355e-03],
                          [3.524609483702912711e-05, -2.033979252945443178e-04, 9.407184674376763221e-05,
                           6.420191853683136557e-04, 9.708296336309190505e-03],
                          [5.448554694265199405e-05, -6.101548506709004874e-04, 2.350679290365150376e-03,
                           -3.885635011674578101e-03, 1.159583515524542641e-02],
                          [2.509069358741303755e-04, -2.745833793134813943e-03, 1.058203324849244474e-02,
                           -1.699644834767223903e-02, 1.990212497011066908e-02],
                          [-3.704820516519310788e-06, 7.131632040647030903e-07, 4.317558143171664072e-04,
                           -2.032816000365457454e-03, 1.143980445336633886e-02],
                          [-4.204569181146790508e-04, 5.206697969508867867e-03, -2.235770446643961518e-02,
                           3.812459558982908558e-02, -9.550759026775079361e-03],
                          [2.231387858012323557e-04, -2.676319546734843950e-03, 1.109950108329610212e-02,
                           -1.855240578227829173e-02, 1.980410645236000550e-02],
                          [-6.691732488527129706e-04, 8.524139067866563085e-03, -3.781911863658569267e-02,
                           6.672331578414059106e-02, -2.545667183248455714e-02]]))

floor_noise_polynomial_coeffs = list(map(numpy.poly1d,[[-5.048115371963505290e-04, 6.443487332502418306e-03, -2.755822559258705451e-02, 4.678398234393036509e-02, -1.495785089947668928e-02],
                     [-9.545526834978516374e-04, 1.093581811555647404e-02, -4.373225745838945494e-02, 7.253724642323591820e-02, -2.947871228053306966e-02],
                     [-5.455072723927560088e-04, 6.747101463631922924e-03, -2.909210838740067598e-02, 5.142524170577621179e-02, -1.844382811363584618e-02],
                     [-9.995979644957279312e-04, 1.143469060900442931e-02, -4.569263627204219608e-02, 7.582041541078338165e-02, -3.164204171591439296e-02],
                     [-8.863856077704932801e-04, 1.035449654357211390e-02, -4.215516396825260553e-02, 7.013816214188291209e-02, -2.774640409510600869e-02],
                     [-7.589651808132531589e-04, 8.863025853651378821e-03, -3.654472453878408050e-02, 6.347248692449157514e-02, -2.594166954960746077e-02],
                     [-8.504010051867004396e-04, 9.925324762449056273e-03, -4.048905515657710147e-02, 6.772132945392388137e-02, -2.613393087775099655e-02],
                     [-8.090098982174978495e-04, 9.212813231610147979e-03, -3.672620889651094511e-02, 6.121778503530301424e-02, -2.345527906678999269e-02],
                     [-7.353309693964817833e-04, 8.579048904274454504e-03, -3.528882623908410898e-02, 6.047480810406839769e-02, -2.345062902295019097e-02],
                     [-1.067629650668973406e-03, 1.218703066328317465e-02, -4.849669060719514341e-02, 8.025101291208461274e-02, -3.434111000309520573e-02],
                     [-6.620368710142958066e-04, 7.804452173703646869e-03, -3.256030413264018875e-02, 5.701722179959860942e-02, -2.212834531433303942e-02],
                     [-8.027046984013080358e-04, 9.059668589149815626e-03, -3.577650957903925205e-02, 5.977638796857074471e-02, -2.403865570972986074e-02],
                     [-7.756165064623434403e-04, 9.040058342914330733e-03, -3.709471264561637477e-02, 6.351020095704701385e-02, -2.495784213739906179e-02],
                     [-9.148838672119842053e-04, 1.053486479175496163e-02, -4.244120168212644345e-02, 7.176011413054944610e-02, -3.031679372554601029e-02],
                     [-1.116809761353775021e-03, 1.287132094085439125e-02, -5.173227197391310622e-02, 8.565007483821840406e-02, -3.605725799412751176e-02],
                     [-7.703254235614598669e-04, 8.875516605196056424e-03, -3.593367461161722787e-02, 6.230828575988874019e-02, -2.710633697060860214e-02]]))


def sample_noise(distance, laser_id, noise_mu=0.0):
    laser_verticals = laser_angles_vlp16
    laser_ids = range(len(laser_angles_vlp16))

    poly = wall_noise_polynomials[laser_id]
    # clamp the distance for now as the polynomial
    # outside of the range measured isn't useful
    if distance < 1.0:
        distance = 1.0
    elif distance > 5.0:
        distance = 5.0

    sampled_sigma = poly(distance)

    return random.gauss(noise_mu, sampled_sigma)


# If the laser noise has to be truely randomize, call this function prior
# to every scan
def randomize_distance_bias(scanner_object, noise_mu = 0.0, noise_sigma = 0.04):
    if scanner_object.velodyne_noise_type == "gaussian":
      for idx in range(len(laser_noise)):
        laser_noise[idx] = random.gauss(noise_mu, noise_sigma)
    elif scanner_object.velodyne_noise_type == "laplace":
      for idx in range(len(laser_noise)):
        laser_noise[idx] = numpy.random.laplace(noise_mu, noise_sigma)
    else:
      raise ValueError("Noise type not supported")

"""
@param world_transformation The transformation for the resulting pointcloud

"""
def scan_advanced(scanner_object,
                  rotation_speed = 10.0,
                  simulation_fps=24,
                  angle_resolution = 0.1728,
                  max_distance = 120,
                  evd_file=None,
                  noise_mu=0.0,
                  noise_sigma=0.03,
                  start_angle = 0.0,
                  end_angle = 360.0,
                  evd_last_scan=True,
                  add_blender_mesh = False,
                  add_noisy_blender_mesh = False,
                  frame_time = (1.0 / 24.0),
                  simulation_time = 0.0,
                  world_transformation=Matrix(),
                  output_laser_id_as_color=False,
                  apply_vertical_correction=False,
                  add_beam_divergence=False):
    
    # First, get the angles for the lasers, depends on the sensor model
    scanner_angles = laser_angles
    # the scanner noise variable isn't used anywhere...
    scanner_noise = laser_noise

    if scanner_object.velodyne_model == BLENSOR_VELODYNE_HDL32E:
      scanner_angles = laser_angles_32

    if scanner_object.velodyne_model == BLENSOR_VELODYNE_VLP16:
      scanner_angles = laser_angles_vlp16

    # These flags control whether to invert the X, Y or Z coordinate
    inv_scan_x = scanner_object.inv_scan_x
    inv_scan_y = scanner_object.inv_scan_y
    inv_scan_z = scanner_object.inv_scan_z    

    start_time = time.time()

    current_time = simulation_time
    delta_rot = angle_resolution*math.pi/180

    evd_storage = evd.evd_file(evd_file)

    xaxis = Vector([1,0,0])
    yaxis = Vector([0,1,0])
    zaxis = Vector([0,0,1])

    rays = []
    ray_info = []

    steps_per_rotation = 360.0/angle_resolution
    time_per_step = (1.0 / rotation_speed) / steps_per_rotation
    angles = end_angle-start_angle
  
    lines = (end_angle-start_angle)/angle_resolution

    horizontal_beam_divergence = 0.18 # about 3.0mrad (0.18deg)
    vertical_beam_divergence = 0.09 # about 1.5 mrad (0.09deg)

    if not add_beam_divergence:
        horizontal_beam_divergence = 0.0
        vertical_beam_divergence = 0.0

    # make a ray for each laser
    ray = Vector([0.0,0.0,0.0])
    for line in range(int(lines)):
        for laser_idx in range(len(scanner_angles)):
            ray.xyz = [0,0,max_distance]

            rot_angle = 1e-6 + start_angle+float(line)*angle_resolution + 180.0
            timestamp = ( (rot_angle-180.0)/angle_resolution) * time_per_step 
            rot_angle = rot_angle%360.0
            ray_info.append([deg2rad(rot_angle), deg2rad(scanner_angles[laser_idx]), timestamp])

            x_rotation = deg2rad(-scanner_angles[laser_idx] + random.uniform(-(vertical_beam_divergence/2), (vertical_beam_divergence/2)))
            y_rotation = deg2rad(rot_angle + random.uniform(-(horizontal_beam_divergence/2), (horizontal_beam_divergence/2)))
            z_rotation = 0.0
            
            rotator = Euler([x_rotation, y_rotation, z_rotation])
            ray.rotate( rotator )
            rays.extend([ray[0],ray[1],ray[2]])

    returns = blensor.scan_interface.scan_rays(rays, max_distance, inv_scan_x = inv_scan_x, inv_scan_y = inv_scan_y, inv_scan_z = inv_scan_z)

    # not sure what this was:
#    for idx in range((len(rays)//3)):
    
    reusable_4dvector = Vector([0.0,0.0,0.0,0.0])
    # structure of one return
    # [<distance>, <x-coord>, <y-coord>, <z-coord>, <object-id>, <rgb-material>, <index>]
    
    for i in range(len(returns)):
        idx = returns[i][-1]
        laser_id = idx % len(scanner_angles)
        x = returns[i][1]
        y = returns[i][2]
        z = returns[i][3]

        if apply_vertical_correction:
            z = z + (vertical_corrections_vlp16[laser_id] / 1000)

        reusable_4dvector.xyzw = (x, y, z, 1.0)

        # this is used for the non-noisy xyz coordinates of point
        vt = (world_transformation * reusable_4dvector).xyz

        # this is used to get the noisy xyz coordinates of point
        v = [x, y, z]

        # make distance noise based on the laser noise array and some gaussian
        distance_noise = laser_noise[laser_id] + random.gauss(noise_mu, noise_sigma)

        # vector length obviously represents the distance between the scanned point and the sensor
        # note: this is already present at returns[i][0] so why is it computed again??
        vector_length = math.sqrt(v[0]**2+v[1]**2+v[2]**2)

        # normalize the vector
        norm_vector = [v[0]/vector_length, v[1]/vector_length, v[2]/vector_length]

        # add the distance noise to the vector length
        vector_length_noise = vector_length+distance_noise

        # multiply each component of the normalized vector with the same distance noise
        reusable_4dvector.xyzw=[norm_vector[0]*vector_length_noise, norm_vector[1]*vector_length_noise, norm_vector[2]*vector_length_noise,1.0]

        # this is the same transformed vector as "vt", except with the noise added
        v_noise = (world_transformation * reusable_4dvector).xyz

        if output_laser_id_as_color:
            color = (laser_id / len(scanner_angles), laser_id / len(scanner_angles), laser_id / len(scanner_angles))
        else:
            color = returns[i][5]

        evd_storage.addEntry(timestamp = ray_info[idx][2], yaw =(ray_info[idx][0]+math.pi)%(2*math.pi), pitch=ray_info[idx][1], distance=vector_length, distance_noise=vector_length_noise, x=vt[0], y=vt[1], z=vt[2], x_noise=v_noise[0], y_noise=v_noise[1], z_noise=v_noise[2], object_id=returns[i][4], color=color)


    current_angle = start_angle+float(float(int(lines))*angle_resolution)

    pre_write_time = time.time()
            
    if evd_file:
        evd_storage.appendEvdFile()

    if not evd_storage.isEmpty():
        scan_data = numpy.array(evd_storage.buffer)
        additional_data = None
        if scanner_object.store_data_in_mesh:
            additional_data = evd_storage.buffer

        if add_blender_mesh:
            # indices 5, 6, 7 are where the x,y,z coordinates are stored
            mesh_utils.add_mesh_from_points_tf(scan_data[:,5:8], "Scan", world_transformation, buffer=additional_data)

        if add_noisy_blender_mesh:
            # indices 8, 9, 10 are where the noisy x,y,z coordinates are stored
            mesh_utils.add_mesh_from_points_tf(scan_data[:,8:11], "NoisyScan", world_transformation, buffer=additional_data) 
            
        bpy.context.scene.update()

    end_time = time.time()
    scan_time = pre_write_time-start_time
    total_time = end_time-start_time
    print ("Elapsed time: %.3f (scan: %.3f)"%(total_time, scan_time))

    return True, current_angle, scan_time




# This Function creates scans over a range of frames

def scan_range(scanner_object,
               frame_start,
               frame_end,
               filename="/tmp/landscape.evd",
               frame_time = (1.0/24.0),
               rotation_speed = 10.0,
               add_blender_mesh=False,
               add_noisy_blender_mesh=False,
               angle_resolution = 0.1728,
               max_distance = 120.0,
               noise_mu = 0.0,
               noise_sigma= 0.02,
               last_frame = True,
               world_transformation=Matrix(),
               output_laser_id_as_color=False,
               apply_vertical_correction=False,
               add_beam_divergence=False):
    start_time = time.time()

    angle_per_second = 360.0 * rotation_speed
    angle_per_frame = angle_per_second * frame_time
    
    try:
        for i in range(frame_start,frame_end):
                bpy.context.scene.frame_current = i

                # Randomize the noise levels on every scan
                randomize_distance_bias(scanner_object)

                ok,start_radians,scan_time = scan_advanced(scanner_object, 
                    rotation_speed=rotation_speed, angle_resolution = angle_resolution, 
                    start_angle = float(i)*angle_per_frame, 
                    end_angle=float(i+1)*angle_per_frame, evd_file = filename, 
                    evd_last_scan=False, add_blender_mesh=add_blender_mesh, 
                    add_noisy_blender_mesh=add_noisy_blender_mesh, 
                    frame_time=frame_time, simulation_time = float(i)*frame_time,
                    max_distance=max_distance, noise_mu = noise_mu, 
                    noise_sigma=noise_sigma, world_transformation=world_transformation,
                    output_laser_id_as_color=output_laser_id_as_color,
                    apply_vertical_correction=apply_vertical_correction,
                    add_beam_divergence=add_beam_divergence)

                if not ok:
                    break
    except BaseException as e:
        print(f"Scan aborted. {e.with_traceback()}")


    if last_frame:
        #TODO: move this into the evd module
        evd_file = open(filename,"a")
        evd_file.buffer.write(struct.pack("i",-1))
        evd_file.close()

    with open(f"{os.path.dirname(filename)}/info.txt", "w") as f:
        info = f"Virtual blensor scan" \
               f"Time of scan: {time.time()}" \
               f"Frame start: {frame_start}" \
               f"Frame end: {frame_end}" \
               f"Frame time: {frame_time}" \
               f"Rotation speed: {rotation_speed}" \
               f"Angle resolution: {angle_resolution}" \
               f"Max distance: {max_distance}" \
               f"Noise mu: {noise_mu}" \
               f"Noise sigma: {noise_sigma}" \
               f"Output laser id as color: {output_laser_id_as_color}" \
               f"Apply vertical correction: {apply_vertical_correction}" \
               f"Add beam divergence: {add_beam_divergence}"
        f.write(info)


    end_time = time.time()
    print ("Total scan time: %.2f"%(end_time-start_time))





