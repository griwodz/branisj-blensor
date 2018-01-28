import math
import sys
import os
import struct 
import ctypes
import time
import random
import bpy
import bmesh
from mathutils import Vector, Euler, Matrix

from blensor import evd
import blensor


def tuples_to_list(tuples):
    l = []
    for t in tuples:
        l.extend(t)
    return l


# Creates and adds a mesh_object from a set of points or from a flattened 
# list
def add_mesh_from_points(points, name="mesh"):
    flattened = points

    if len(points) > 0:
        if type(points[0]) == type([]) or type(points[0]) == type(()):
            flattened = tuples_to_list(points)            

    mesh = bpy.data.meshes.new(name+"_mesh")
    mesh.vertices.add(len(flattened)//3)
    mesh.vertices.foreach_set("co",flattened)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.scene.objects.link(obj)

    scanner = bpy.context.object
    if scanner.show_in_frame: 
      blensor.show_in_frame(obj, bpy.context.scene.frame_current)



"""Add a mesh from points or from a flattened list and transform it according to
   world_transformation, or according to the transformation of the camera
"""
def add_mesh_from_points_tf(points, name="Scan", world_transformation = Matrix(), buffer = None):
    bm = bmesh.new()
    for p in points:
        bm.verts.new(p)

    bm.verts.ensure_lookup_table()

    if buffer:
        color_red_id = bm.verts.layers.float.new("color_red")
        color_green_id = bm.verts.layers.float.new("color_green")
        color_blue_id = bm.verts.layers.float.new("color_blue")
        for idx,c in enumerate(buffer):
            v = bm.verts[idx]
            v[color_red_id] = c[12]
            v[color_green_id] = c[13]
            v[color_blue_id] = c[14]

    mesh = bpy.data.meshes.new(name+"_mesh")
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    mesh_object = bpy.data.objects.new("{0}.{1}".format(name,bpy.context.scene.frame_current), mesh)
    bpy.context.scene.objects.link(mesh_object)

    scanner = bpy.context.object
    if scanner.show_in_frame: 
      blensor.show_in_frame(mesh_object, bpy.context.scene.frame_current)

    if world_transformation == Matrix():
       mesh_object.matrix_world = bpy.context.object.matrix_world



