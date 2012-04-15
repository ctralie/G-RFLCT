from Primitives3D import *
from Graphics3D import *
from Shapes3D import *
from PolyMesh import *
from EMScene import *
from Cameras3D import *
import math
import Image

def ConstructRayThroughPixel(camera, x, y, width, height):
	tanY = math.tan(camera.yfov)
	tanX = math.tan(camera.yfov*float(width)/float(height))
	P0 = camera.eye
	towards = camera.towards
	up = camera.up
	right = towards % up
	towards.normalize()
	up.normalize()
	right.normalize()
	P1 = P0 + towards - tanX*right
	P2 = P0 + towards + tanX*right
	P = P1 + ((float(x) + 0.5) / float(width))*(P2 - P1)
	farwidth = (P2 - P1).Length()
	P = P + (float(y-float(height)/2.0)/(float(height)/2.0)*tanY)*up
	directionVec = P - P0
	directionVec.normalize()
	return Ray3D(P0, directionVec)
	
def RayTraceImage(scene, camera, width, height, filename):
	im = Image.new("RGB", (width, height))
	pix = im.load()
	for x in range(0, width):
		print "Finished Row %i"%x
		for y in range(0, height):
			ray = ConstructRayThroughPixel(camera, x, y, width, height)
			intersection = scene.getRayIntersection(ray)
			if intersection != None:
				pix[x, y] = (255, 255, 255)
			else:
				pix[x, y] = (20, 20, 20)
	im.save(filename)
