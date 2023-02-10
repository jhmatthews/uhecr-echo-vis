import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

def run():
	'''
	simple script to combine the movie frames into one summary frame
	'''
	frame_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Movies/Frames/'))
	Models = ["ModelA","ModelB","ModelC"]
	for Model in Models:
		time_stamps = np.arange(0,87,1)
		print ("summary plot for {}".format(Model))
		for t in tqdm(time_stamps):

			plt.figure()
			xy_fname = "{}/xy_{}_{:03d}.png".format(frame_dir, Model, t)
			map_fname = "{}/map_{}_{:d}_sm20.png".format(frame_dir, Model, t)
			time_fname = "{}/timeseries_{}_{:03d}.png".format(frame_dir, Model, t)

			plt.subplot(2,2,1)
			im = plt.imread(xy_fname)
			plt.imshow(im)
			plt.axis("off")

			plt.subplot(2,2,2)
			im = plt.imread(time_fname)
			plt.imshow(im)
			plt.axis("off")

			plt.subplot(2,1,2)
			im = plt.imread(map_fname)
			plt.imshow(im)
			plt.axis("off")

			plt.subplots_adjust(hspace=0, wspace=0,right=1,left=0,bottom=0,top=0.98)
			plt.savefig("{}/summary_{}_{:03d}.png".format(frame_dir, Model, t), dpi=600)
			plt.close("all")
		
if __name__ == "__main__":
	run()
