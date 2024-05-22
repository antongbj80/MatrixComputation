import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
 
def pic_compress(r,pic_array):
    global u,sigma,vt,sig,new_pic
 
    u,sigma,vt = np.linalg.svd(pic_array)
    sig = np.diag(sigma[: r])
    new_pic = np.dot(np.dot(u[:, :r],sig),vt[:r,:])
    size = u.shape[0] * r +sig.shape[0] * sig.shape[1]+ r * vt.shape[1]
    print("compress size r = {}: {}".format(r, size))
    return new_pic,size
 
img = Image.open('lecture 4/1gray.png')
ori_img = np.array(img)
print("original size:"+str(ori_img.shape[0]*ori_img.shape[1]))

new_img_5, size  = pic_compress(5,ori_img)
new_img_10, size  = pic_compress(10,ori_img)
new_img_20, size  = pic_compress(20,ori_img)
new_img_30, size  = pic_compress(30,ori_img)
new_img_50, size  = pic_compress(50,ori_img)

fig, ax = plt.subplots(2, 3, figsize=(15, 10))
ax[0, 0].imshow(ori_img, cmap='gray')
ax[0, 0].set_title("before compress")

ax[0, 1].imshow(new_img_5, cmap='gray')
ax[0, 1].set_title("after compress_r = 5")

ax[0, 2].imshow(new_img_10, cmap='gray')
ax[0, 2].set_title("after compress_r = 10")

ax[1, 0].imshow(new_img_20, cmap='gray')
ax[1, 0].set_title("after compress_r = 20")

ax[1, 1].imshow(new_img_30, cmap='gray')
ax[1, 1].set_title("after compress_r = 30")

ax[1, 2].imshow(new_img_50, cmap='gray')
ax[1, 2].set_title("after compress_r = 50")

plt.tight_layout(pad=3.0, w_pad=0.5, h_pad=3.0)
plt.show()