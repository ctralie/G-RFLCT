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
	P1 = towards - tanX*right
	P2 = towards + tanX*right
	directionVec = P1 + (float(x)/float(width))*(P2 - P1)
	farwidth = (P2 - P1).Length()
	directionVec = directionVec + (float(height-y)/float(height)-0.5)*(tanY*up)
	directionVec.normalize()
	return Ray3D(P0, directionVec)
	
def RayTraceImage(scene, camera, width, height, filename):
	rayPoints = []
	rayNormals = []
	im = Image.new("RGB", (width, height))
	pix = im.load()
	for x in range(0, width):
		print "Finished Row %i"%x
		for y in range(0, height):
			ray = ConstructRayThroughPixel(camera, x, y, width, height)
			intersection = scene.getRayIntersection(ray)
			if intersection != None:
				pix[x, y] = (255, 255, 255)
				rayPoints.append(ray.P0)
				rayPoints.append(intersection[1])
				rayNormals.append(intersection[1])
				rayNormals.append(intersection[1]+0.1*intersection[2])
				print "Face %i"%intersection[3].ID
			else:
				pix[x, y] = (20, 20, 20)
	im.save(filename)
	return (rayPoints, rayNormals)

if __name__ == '__main__':
	im = Image.new("RGB", (400, 400))
	pix = im.load()
	for i in range(0, 400):
		pix[i, i] = (255, 255, 255)
	im.save("diag.png")
