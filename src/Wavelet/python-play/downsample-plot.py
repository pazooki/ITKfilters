import numpy as np
# import matplotlib
import matplotlib.pyplot as plt

xspacing = 0.04
x = np.arange(0., 6 * np.pi, xspacing)
print x.size
plt.subplot(3, 2, 1)
# y = np.sin(x/2.)
y = np.sin(x/2.) + np.sin(2.*x) + np.sin(16. * x)
plt.plot(x, y, 'b.')
plt.subplot(3, 2, 2)
freq = np.fft.fftfreq(x.shape[-1])
fft = np.fft.fft(y)
plt.plot(freq, fft)
# DOWNSAMPLED
plt.subplot(3, 2, 3)
x2 = x[::2].copy()
y2 = np.sin(x2/2.) + np.sin(2.*x2) + np.sin(16. * x2)
print x2.size
plt.plot(x2, y2, 'r+')
plt.subplot(3, 2, 4)
freq2 = np.fft.fftfreq(x2.shape[-1])
fft2 = np.fft.fft(y2)
plt.plot(freq2, fft2)

# UPSAMPLED
xspacingu = 0.02
xu = np.arange(0., 6 * np.pi, xspacingu)
print xu.size
plt.subplot(3, 2, 5)
yu = np.sin(xu/2.) + np.sin(2.*xu) + np.sin(16. * xu)
plt.plot(xu, yu, 'r+')
plt.subplot(3, 2, 6)
frequ = np.fft.fftfreq(xu.shape[-1])
fftu = np.fft.fft(yu)
plt.plot(frequ, fftu)


plt.show()
