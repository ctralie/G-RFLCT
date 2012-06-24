from Primitives3D import *
from Graphics3D import *
from Shapes3D import *
from PolyMesh import *
from EMScene import *
from Cameras3D import *
import math
import Image

def ConstructRayThroughPixel(camera, xScale, yScale, x, y, width, height):
	P0 = camera.eye
	towards = camera.towards
	up = camera.up
	right = towards % up
	towards.normalize()
	up.normalize()
	right.normalize()
	xCoord = float(x)/float(width) - 0.5
	yCoord = float(y)/float(height) - 0.5
	directionVec = towards + xCoord*2*xScale*right + yCoord*2*yScale*up
	directionVec.normalize()
	return Ray3D(P0, directionVec)
	
def RayTraceImage(scene, camera, width, height, filename):
	yfov = camera.yfov
	xfov = yfov*float(width)/float(height)
	xScale = math.tan(xfov/2)
	yScale = math.tan(yfov/2)
	print "xScale = %g, yScale = %g\n"%(xScale, yScale)
	rayPoints = []
	rayNormals = []
	im = Image.new("RGB", (width, height))
	pix = im.load()
	for x in range(0, width):
		print "Finished Row %i"%x
		for y in range(0, height):
			ray = ConstructRayThroughPixel(camera, xScale, yScale, x, y, width, height)
			intersection = scene.getRayIntersection(ray)
			flipY = height - y - 1
			if intersection != None:
				pix[x, flipY] = (255, 255, 255)
				rayPoints.append(ray.P0)
				rayPoints.append(intersection[1])
				rayNormals.append(intersection[1])
				rayNormals.append(intersection[1]+0.1*intersection[2])
			else:
				pix[x, flipY] = (20, 20, 20)
	im.save(filename)
	return (rayPoints, rayNormals)

if __name__ == '__main__':
	im = Image.new("RGB", (400, 400))
	pix = im.load()
	for i in range(0, 400):
		pix[i, i] = (255, 255, 255)
	im.save("diag.png")
