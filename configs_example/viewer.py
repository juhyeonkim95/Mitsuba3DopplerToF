import matplotlib.pyplot as plt
import cv2

image = cv2.imread("scene.exr",  cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
image = image * 0.0015
plt.imshow(image[:, :, 0], vmin=-1e-6, vmax=1e-6)
plt.show()