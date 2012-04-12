class RGB3D(object):
	def __init__(self, R = 1.0, G = 1.0, B = 1.0):
		[self.R, self.G, self.B] = [R, G, B]

class OpticalMaterial(object):
	def __init__(self, ka = RGB3D(0.2, 0.2, 0.2), kd = RGB3D(0.5, 0.5, 0.5), ks = RGB3D(0.5, 0.5, 0.5), kt = RGB3D(0.0, 0.0, 0.0), emission = RGB3D(0.0, 0.0, 0.0)):
		[self.ka, self.kd, self.ks] = [ka, kd, ks]
		[self.kt, self.emission] = [kt, emission]
